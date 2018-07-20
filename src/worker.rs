use std::collections::VecDeque;
use std::sync::{Arc, Mutex, Condvar};
use std::thread;
use std::thread::JoinHandle;

#[derive(Debug, new)]
struct Worker {
    index: usize,
    thread: JoinHandle<()>,
}

pub struct WorkerPool<D: Send + 'static> {
    workers: Vec<Worker>,
    work_queue: Arc<(Mutex<VecDeque<D>>, Condvar)>,
    num_waiting_workers: Arc<(Mutex<usize>, Condvar)>,
}

impl<D: Send> WorkerPool<D> {
    pub fn new(num_workers: usize, processor: Box<Fn(D) + Send + Sync>) -> Self {
        let work_queue = Arc::new((Mutex::new(VecDeque::new()), Condvar::new()));
        let processor = Arc::new(processor);
        let num_waiting_workers = Arc::new((Mutex::new(0), Condvar::new()));

        let workers = (0..num_workers).map(|i| {
            let work_queue2 = work_queue.clone();
            let processor2 = processor.clone();
            let num_waiting_workers2 = num_waiting_workers.clone();

            let thread = thread::spawn(move || {
                loop {
                    // Here, it is assumed that the thread will have to wait on the condvar of
                    // the queue. Therefore the number of waiting workers is incremented. Since
                    // the owning thread might wait on the condvar of the num_waiting_workers 
                    // variable, it has to be notified.
                    {
                        let (ref num_waiting_workers, ref condvar) = *num_waiting_workers2;
                        let mut num_waiting_workers = num_waiting_workers.lock().unwrap(); // @TODO: Handle the unwrap
                        *num_waiting_workers += 1;

                        condvar.notify_one();
                    }

                    {
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
                            queue = condvar.wait(queue).unwrap();
                        }
                    }

                    // After the worker gets new work and starts processing it, the number of
                    // waiting workers has to be decremented.
                    {
                        let (ref num_waiting_workers, ref condvar) = *num_waiting_workers2;
                        let mut num_waiting_workers = num_waiting_workers.lock().unwrap(); // @TODO: Handle the unwrap
                        *num_waiting_workers -= 1;

                        condvar.notify_one();
                    }

                    {
                        // @TODO: Reuse the queue variable from above and test whether this lead
                        // to the thread being blocked too long.
                        let (ref queue, ..) = *work_queue2;
                        let mut queue = queue.lock().unwrap(); // @TODO: Handle the unwrap
                        if let Some(work) = queue.pop_front() {
                            // Not dropping the guard has the consequence to blocking parallelism altogether.
                            drop(queue);
                            processor2(work);
                        }
                    }
                }
            });

            Worker::new(i, thread)
        }).collect();

        WorkerPool {
            workers: workers,
            work_queue: work_queue,
            num_waiting_workers: num_waiting_workers,
        }
    }

    pub fn queue_len(&self) -> usize {
        let (ref queue, ..) = *self.work_queue;
        let queue = queue.lock().unwrap(); // @TODO: Handle the unwrap
        queue.len()
    }

    pub fn process(&self, job_data: D) {
        let (ref queue, ref condvar) = *self.work_queue;
        let mut queue = queue.lock().unwrap(); // @TODO: Handle the unwrap
        queue.push_back(job_data);

        condvar.notify_all(); // @TODO: Only one should be notified.
    }

    pub fn wait(&self) {
        while {
            // The num_waiting_workers variable is incremented before the worker checks 
            // whether there is work in the queue. When there is no work, the worker
            // waits on the condvar of the queue. In this case the increment does make
            // sense and represents the state of the worker. When there is work in the
            // queue, the worker can go on to process the work. In this case, the variable
            // num_waiting_workers is decremented again. Therefore it is possible (only
            // for a split second) that all workers are considered waiting, although 
            // there is work to do.
            // For this reason, it is necessary to wait on the condvar of num_waiting_workers
            // and check if all workers are waiting BUT ALSO whether the queue is really 
            // empty.
            // Keep in mind that while this function is run, no work can be inserted into
            // the queue as only the owning thread can run the process function. (But not
            // while this one is executed.)

            {
                let (ref num_waiting_workers, ref cvar) = *self.num_waiting_workers;
                let mut num_waiting_workers = num_waiting_workers.lock().unwrap(); // @TODO: Handle the unwrap
                while *num_waiting_workers != self.workers.len() {
                    num_waiting_workers = cvar.wait(num_waiting_workers).unwrap(); // @TODO: Handle the unwrap
                }
            }

            self.queue_len() > 0
        }{}
    }

    fn panic_mutex_lock_message(mutex_string: &str) -> String {
        format!("Could not aquire the mutex for the {} within the WorkerPool as the mutex is poisened. This indicates that some erroneous condition on another thread was not handled correctly.", mutex_string)
    }
}
