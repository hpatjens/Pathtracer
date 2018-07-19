use std::collections::VecDeque;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, Condvar};
use std::thread;
use std::thread::JoinHandle;

#[derive(Clone, Copy, Debug, PartialEq)]
enum WorkerState {
    Working,
    Waiting,
}

#[derive(Debug)]
struct Worker {
    index: usize,
    state: Arc<Mutex<WorkerState>>,
    thread: JoinHandle<()>,
}

impl Worker {
    pub fn new(index: usize, state: Arc<Mutex<WorkerState>>, thread: JoinHandle<()>) -> Self {
        Worker {
            index: index,
            state: state,
            thread: thread,
        }
    }
}

pub struct WorkerPool<D: Send + 'static> {
    workers: Vec<Worker>,
    work_queue: Arc<(Mutex<VecDeque<D>>, Condvar)>,
    processor: Arc<Mutex<Box<FnMut(D) + Send>>>,
    stop_request: AtomicBool,
}

impl<D: Send> WorkerPool<D> {
    pub fn new(num_workers: usize, processor: Box<FnMut(D) + Send>) -> Self {
        let work_queue = Arc::new((Mutex::new(VecDeque::new()), Condvar::new()));
        let processor = Arc::new(Mutex::new(processor));

        let workers = (0..num_workers).map(|i| {
            let work_queue2 = work_queue.clone();
            let processor2 = processor.clone();

            let worker_state = Arc::new(Mutex::new(WorkerState::Waiting));
            let worker_state2 = worker_state.clone();

            let thread = thread::spawn(move || {
                loop {
                    let (ref queue, ref condvar) = *work_queue2;
                    let mut queue = match queue.lock() {
                        Ok(queue) => queue,
                        Err(_poisened_guard) => {
                            // When another thread panics while having locked the mutex. It is
                            // considered poisened and cannot be aquired. However, Err returns the
                            // poisened guard which could be used anyway which is not a good idea 
                            // in this case as the state of the queue is unknown.
                            panic!(Self::panic_mutex_lock_message("work_queue"));
                        },
                    };
                    while queue.len() <= 0 {
                        {
                            // @TODO: Make a function that encapsulates this.
                            // @TODO: Handle the unwrap()
                            let mut worker_state = worker_state2.lock().unwrap();
                            *worker_state = WorkerState::Waiting;
                        }
                        queue = condvar.wait(queue).unwrap();
                    }
                    {
                        // @TODO: Make a function that encapsulates this.
                        // @TODO: Handle the unwrap()
                        let mut worker_state = worker_state2.lock().unwrap();
                        *worker_state = WorkerState::Waiting;
                    }
                    
                    if let Some(work) = queue.pop_front() {
                        let mut processor = match processor2.lock() {
                            Ok(processor) => processor,
                            Err(_poisened_guard) => panic!(Self::panic_mutex_lock_message("processor")),
                        };
                        (&mut *processor)(work);
                    }
                }
            });

            Worker::new(i, worker_state, thread)
        }).collect();

        WorkerPool {
            workers: workers,
            work_queue: work_queue,
            processor: processor,
            stop_request: AtomicBool::new(false),
        }
    }

    pub fn set_processor(&self, new_processor: Box<FnMut(D) + Send>) {
        let mut processor = self.processor.lock().unwrap();
        *processor = new_processor;
    }

    pub fn queue_len(&self) -> usize {
        let (ref queue, ..) = *self.work_queue;
        let queue = queue.lock().unwrap();
        queue.len()
    }

    pub fn process(&self, job_data: D) -> bool {
        if self.stop_request.load(Ordering::Relaxed) {
            return false;
        }

        let (ref queue, ref condvar) = *self.work_queue;
        let mut queue = match queue.lock() {
            Ok(queue) => queue,
            Err(_poisened_guard) => panic!(Self::panic_mutex_lock_message("work_queue")),
        };
        queue.push_back(job_data);

        condvar.notify_all();

        true
    }

    pub fn wait(&self) {
        self.stop_request.store(true, Ordering::Relaxed);

        // @TODO: Replace the busy waiting with a condvar
        loop {
            let queue_empty = {
                let (ref queue, ..) = *self.work_queue;
                let queue = queue.lock().unwrap();
                queue.len() == 0
            };
            if queue_empty && self.all_workers_waiting() {
                break;
            }
        }

        self.stop_request.store(false, Ordering::Relaxed);
    }

    pub fn all_workers_waiting(&self) -> bool {
        let mut all_waiting = true;
        for worker in &self.workers {
            // @TODO: Handle the unwrap
            let state = worker.state.lock().unwrap();
            if *state == WorkerState::Working {
                all_waiting = false;
            }
        }
        all_waiting
    }

    fn panic_mutex_lock_message(mutex_string: &str) -> String {
        format!("Could not aquire the mutex for the {} within the WorkerPool as the mutex is poisened. This indicates that some erroneous condition on another thread was not handled correctly.", mutex_string)
    }
}
