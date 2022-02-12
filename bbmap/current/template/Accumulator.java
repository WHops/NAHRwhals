package template;

/**
 * Interface for accumulating statistics captured by threads.
 * 
 * @author Brian Bushnell
 * @date November 19, 2015
 *
 * @param <T>
 */
public interface Accumulator<T> {
	
	/** Accumulate personal variables */
	public void accumulate(T t);
	
	/** True if it finished successfully */
	public boolean success();
	
}
