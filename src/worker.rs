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
                        let mut num_waiting_workers = num_waiting_workers.lock().expect("Could not aquire the mutex for incrementing the waiting workers.");
                        *num_waiting_workers += 1;

                        condvar.notify_one();
                    }

                    // At this point the worker is considered waiting although there might be 
                    // work within the queue. Due to this inaccuracy the wait method has to
                    // check whether the queue is really empty when all workers signal that
                    // they are waiting.
                    {
                        let (ref queue, ref condvar) = *work_queue2;
                        let mut queue = queue.lock().expect("Could not aquire the mutex for waiting on the queue.");
                        while queue.len() <= 0 {
                            queue = condvar.wait(queue).expect("Could not aquire the mutex for waiting on the queue while waiting on the condvar.");
                        }
                    }

                    // After the worker gets new work and starts processing it, the number of
                    // waiting workers has to be decremented.
                    {
                        let (ref num_waiting_workers, ref condvar) = *num_waiting_workers2;
                        let mut num_waiting_workers = num_waiting_workers.lock().expect("Could not aquire the mutex for decrementing the waiting workers.");
                        *num_waiting_workers -= 1;

                        condvar.notify_one();
                    }

                    // Here, the lock to the queue is aquired again. It would be possible to
                    // keep the lock from above but this would reduce the chance for other
                    // threads to access the queue.
                    // It is important to drop the lock with in the if let so that the mutex
                    // for the queue is unlocked while processing the work. Otherwise parallelism
                    // would not be possible while processing.
                    {
                        let mut queue = work_queue2.0.lock().expect("Could not aquire the mutex for dequeuing");
                        if let Some(work) = queue.pop_front() {
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
        self.work_queue.0.lock().expect("Could not aquire the mutex to get the queue length.").len()
    }

    pub fn process(&self, job_data: D) {
        let (ref queue, ref condvar) = *self.work_queue;
        let mut queue = queue.lock().expect("Could not aquire the mutex of the queue for inserting.");
        queue.push_back(job_data);

        condvar.notify_one();
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
                let (ref num_waiting_workers, ref condvar) = *self.num_waiting_workers;
                let mut num_waiting_workers = num_waiting_workers.lock().expect("Could not aquire the mutex for num_waiting_workers in wait.");
                while *num_waiting_workers != self.workers.len() {
                    num_waiting_workers = condvar.wait(num_waiting_workers).expect("Could not aquire the mutex for num_waiting_workers in wait while waiting on the condvar.");
                }
            }

            self.queue_len() > 0
        }{}
    }
}
