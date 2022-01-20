package structures;

import java.util.ArrayList;

public final class Heap<T extends Comparable<? super T>> {
	
	@SuppressWarnings("unchecked")
	public Heap(int maxSize, boolean rollover_){
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=(T[])new Comparable[len];
		rollover=rollover_;
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	public boolean add(T t){
		
		assert(size==0 || array[size]!=null);
		assert(rollover || size<CAPACITY);
		
		if(size>=CAPACITY){
			
			if(t.compareTo(array[1])<=0){return false;}
			
			poll(); //Turns into a rolling buffer by removing smallest value.
			
//			{//This is a more efficient alternative to poll() and percDown(), but the result is slightly different.
//				array[1]=t;
//				percUp(1);
//				return true;
//			}
		}
		assert(size<CAPACITY);
		
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=null);
		size++;
		array[size]=t;
		percDown(size);
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return true;
	}
	
	public T peek(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
//		assert(array[1]==queue.peek()) : size+", "+queue.size()+"\n"+
//			array[1]+"\n"+
//			array[2]+" , "+array[3]+"\n"+
//			array[4]+" , "+array[5]+" , "+array[6]+" , "+array[7]+"\n"+
//			queue.peek()+"\n";
		//assert(testForDuplicates());
		return array[1];
	}
	
	public T poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
		T t=array[1];
//		assert(t==queue.poll());
		array[1]=array[size];
		array[size]=null;
		size--;
		if(size>0){percUp(1);}
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return t;
	}
	
	private void percDown(int loc){
		//assert(testForDuplicates());
		assert(loc>0);
		if(loc==1){return;}
		int next=loc/2;
		T a=array[loc];
		T b=array[next];
		assert(a!=b);
		if(a.compareTo(b)<0){
			array[next]=a;
			array[loc]=b;
			percDown(next);
		}
	}
	
	private void percUp(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		int next1=loc*2;
		int next2=next1+1;
		if(next1>size){return;}
		T a=array[loc];
		T b=array[next1];
		T c=array[next2];
		assert(a!=b);
		assert(b!=c);
		assert(b!=null);
		//assert(testForDuplicates());
		if(c==null || b.compareTo(c)<1){
			if(a.compareTo(b)>0){
				array[next1]=a;
				array[loc]=b;
				//assert(testForDuplicates());
				percUp(next1);
			}
		}else{
			if(a.compareTo(c)>0){
				array[next2]=a;
				array[loc]=c;
				//assert(testForDuplicates());
				percUp(next2);
			}
		}
	}
	
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}
	
	public boolean hasRoom(){
		return size<CAPACITY;
	}
	
	public void clear(){
//		queue.clear();
		for(int i=1; i<=size; i++){array[i]=null;}
		size=0;
	}
	
	public int size(){
		return size;
	}
	
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
			}
		}
		return true;
	}
	
	public ArrayList<T> toList(){
		ArrayList<T> list=new ArrayList<T>(size);
		for(int i=0, lim=size; i<lim; i++){
			list.add(poll());
		}
		assert(isEmpty());
		return list;
	}
	
	private final T[] array;
	private final int CAPACITY;
	private final boolean rollover;
	private int size=0;
	
//	private PriorityQueue<T> queue;
	
}
