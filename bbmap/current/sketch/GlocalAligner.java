package sketch;

import shared.Tools;

public final class GlocalAligner {

	GlocalAligner(){}
	
	public float alignForward(final byte[] query, final byte[] ref){
		float a=alignForwardInner(query, ref);
		float b=alignForwardInner(ref, query);
		return Tools.max(a, b);
	}
	
	/**
	 * @param query
	 * @param ref
	 * @return Identity
	 */
	public float alignForwardInner(final byte[] query, final byte[] ref){
//		if(ref.length<query.length){return alignForward(ref, query);}
		final int big=ref.length;
		final int small=query.length;
//		assert(big>=small);
		final int arrayLength=big;

		//Stack allocated; faster.
		final short[] array1=new short[arrayLength+1], array2=new short[arrayLength+1];
		final short[] consumed1=new short[arrayLength+1], consumed2=new short[arrayLength+1];

		short[] prev=array1, next=array2;
		short[] prevConsumed=consumed1, nextConsumed=consumed2;
		
//		for(int i=0; i<=arrayLength-query.length; i++){prev[i]=0;}
//		for(short i=(short)(arrayLength-query.length), score=0; i<=arrayLength; i++, score+=pointsDel) {
//			prev[i]=score;
//		}
		final int rmax=ref.length-1;
		final int qmax=query.length-1;
		
		int maxScore=0;
		int maxConsumed=0;
		for(short qpos=0; qpos<query.length; qpos++){
//			prev[0]=(short)(pointsIns*qpos);
			
			final byte q=query[qpos];
			for(int rpos=0, apos=1; rpos<ref.length; rpos++, apos++){
				final byte r=ref[rpos];
				final boolean match=(q==r && q!='N');
				final boolean rEnd=(rpos<1 || rpos>=rmax);
				final boolean qEnd=(qpos<1 || qpos>=qmax);
				final short vScore=(short) (prev[apos]+(rEnd ? 0 : pointsIns));
				final short hScore=(short) (next[apos-1]+(qEnd ? 0 : pointsDel));
				final short dScore=(short) ((match ? pointsMatch : pointsSub)+prev[apos-1]);
				
				short score, consumed;
				if(dScore>=vScore && dScore>=hScore){
					score=dScore;
					consumed=(short)(prevConsumed[apos-1]+1);
				}else if(vScore>=hScore){
					score=vScore;
//					nextConsumed[apos]=(short)(prevConsumed[apos]+(rEnd ? 0 : 1));
					consumed=(short)(prevConsumed[apos]);
				}else{
					score=vScore;
					consumed=(short)(nextConsumed[apos-1]);
				}

				if(score<0){
					score=0;
					consumed=0;
				}
				if(score>maxScore){
					maxScore=score;
					maxConsumed=consumed;
				}
				next[apos]=score;
				nextConsumed[apos]=consumed;
//				//Should be branchless conditional moves
//				short score=(dScore>=vScore ? dScore : vScore);
//				score=(hScore>score ? hScore : score);
//				next[apos]=score;
			}
//			iters+=arrayLength;
			
			short[] temp=prev;
			prev=next;
			next=temp;
			temp=prevConsumed;
			prevConsumed=nextConsumed;
			nextConsumed=temp;
		}

//		short maxScore=Short.MIN_VALUE;
//		short maxConsumed=Short.MIN_VALUE;
//		for(short apos=1; apos<prev.length; apos++){//Grab high score from last iteration
//			short score=prev[apos];
//			if(score>=maxScore){
//				maxScore=score;
//				maxConsumed=prevConsumed[apos];
//			}
//		}
		
//		assert(false);
//		System.err.println("maxScore="+maxScore+"; maxConsumed="+maxConsumed);
		if(maxConsumed<400){return 0;}
		
//		int maxPossibleScore=(small*(pointsMatch-pointsSub));
//		int rescaledScore=maxScore-small*pointsSub;
//		final float ratio=rescaledScore/(float)maxPossibleScore;
//		return ratio;
		
		int maxPossibleScore=(maxConsumed*(pointsMatch-pointsSub));
		int rescaledScore=maxScore-maxConsumed*pointsSub;
		final float ratio=rescaledScore/(float)maxPossibleScore;
		return ratio;
	}


	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long iters=0;
	
	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final short pointsMatch = 1;
	public static final short pointsSub = -1;
	public static final short pointsDel = -1;
	public static final short pointsIns = -1;
	
}
