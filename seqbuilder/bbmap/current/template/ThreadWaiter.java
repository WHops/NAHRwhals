package template;

public class ThreadWaiter {
	
	/** Wait for completion of all threads */
	public static final <T extends Thread> boolean waitForThreads(Iterable<T> iter){

		//Wait for completion of all threads
		boolean success=true;
		for(T t : iter){

			//Wait until this thread has terminated
			while(t.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					t.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
		}
		
		return success;
	}
	
	public static final <T extends Thread> void startThreads(Iterable<T> iter){
		for(Thread t : iter){t.start();}
	}
	
	/**
	 * @param iter List of Threads.
	 * @return success
	 */
	public static final <T extends Thread> boolean startAndWait(Iterable<T> iter){
		startThreads(iter);
		return waitForThreads(iter);
	}
	
	/**
	 * @param iter List of Threads.
	 * @return success
	 */
	public static final <T extends Thread> boolean startAndWait(Iterable<T> iter, Accumulator<T> acc){
		startThreads(iter);
//		assert(false);
		return waitForThreads(iter, acc);
	}
	
	/** Wait for completion of all threads, and accumulate results */
	public static final <T extends Thread> boolean waitForThreads(Iterable<T> iter, Accumulator<T> acc){

		waitForThreads(iter);
		accumulate(iter, acc);
		return acc.success();
		
	}
	
	/** Accumulate results from all threads */
	private static final <T> boolean accumulate(Iterable<T> iter, Accumulator<T> acc){

		//Wait for completion of all threads
		for(T t : iter){
//			assert(t.getState()==Thread.State.TERMINATED);//Not strictly necessary; requires T to be a thread.

			//Accumulate per-thread statistics
			acc.accumulate(t);
		}
		
		return acc.success();
	}
	
}
