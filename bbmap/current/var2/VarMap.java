package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

public class VarMap implements Iterable<Var> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/
	
	VarMap(ScafMap scafMap_){
		this(scafMap_, -1, -1, -1, -1, -1);
	}

	@SuppressWarnings("unchecked")
	VarMap(ScafMap scafMap_, int ploidy_, double pairingRate_, double totalQualityAvg_,
			double mapqAvg_, double readLengthAvg_){
		scafMap=scafMap_;
		ploidy=ploidy_;
		properPairRate=pairingRate_;
		totalQualityAvg=totalQualityAvg_;
		totalMapqAvg=mapqAvg_;
		readLengthAvg=readLengthAvg_;
		maps=new ConcurrentHashMap[WAYS];
		for(int i=0; i<WAYS; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
	}
	
//	public static VarMap loadVars(String fname, ScafMap scafMap){
//		final ByteFile bf=ByteFile.makeByteFile(fname, true);
//		final VarMap varMap=new VarMap(scafMap);
//		final byte delimiter='\t';
//		int ploidy=-1;
//		double pairingRate=-1;
//		double mapqAvg=-1;
//		double totalQualityAvg=-1;
//		double readLengthAvg=-1;
//		byte[] line=bf.nextLine();
//		while(line!=null && line.length>0){
//			if(line[0]!='#'){
//				Var v=new Var(line, delimiter);
//				varMap.addUnsynchronized(v);
//			}else{
//				String[] split=new String(line).split("\t");
//				String a=split[0], b=(split.length>1 ? split[1] : null);
//				assert(split.length>1) : new String(line);
//				if(a.equalsIgnoreCase("#ploidy")){
//					ploidy=Integer.parseInt(b);
//				}else if(a.equalsIgnoreCase("#pairingRate")){
//					pairingRate=Double.parseDouble(b);
//				}else if(a.equalsIgnoreCase("#totalQualityAvg")){
//					totalQualityAvg=Double.parseDouble(b);
//				}else if(a.equalsIgnoreCase("#mapqAvg")){
//					mapqAvg=Double.parseDouble(b);
//				}else if(a.equalsIgnoreCase("#readLengthAvg")){
//					readLengthAvg=Double.parseDouble(b);
//				}
//			}
//			line=bf.nextLine();
//		}
//		bf.close();
//		varMap.ploidy=ploidy;
//		varMap.properPairRate=(double)pairingRate;
//		varMap.totalQualityAvg=(double)totalQualityAvg;
//		varMap.totalMapqAvg=(double)mapqAvg;
//		varMap.readLengthAvg=(double)readLengthAvg;
//		return varMap;
//	}
//	
//	//#CHROM POS    ID        REF  ALT     QUAL
//	public static VarMap loadVcf(String fname, ScafMap scafMap){
//		ByteFile bf=ByteFile.makeByteFile(fname, true);
//		VarMap varMap=new VarMap(scafMap);
//		byte[] line=bf.nextLine();
//		while(line!=null && line.length>0){
//			if(line[0]!='#'){
//				Var v;
//				try {
//					v = Var.fromVCF(line, scafMap);
//					varMap.addUnsynchronized(v);
//				} catch (Exception e) {
//					System.err.println("Unable to parse VCF line: '"+new String(line)+"'");
////					throw new RuntimeException(e);
//				}
//			}else{
//				String[] split=new String(line).split("=");
//				if(split.length==2){
//					String a=split[0], b=split[1];
//					if(a.equalsIgnoreCase("##ploidy")){
//						varMap.ploidy=Integer.parseInt(b);
//					}else if(a.equalsIgnoreCase("##properPairRate")){
//						varMap.properPairRate= Double.parseDouble(b);
//					}else if(a.equalsIgnoreCase("##totalQualityAvg")){
//						varMap.totalQualityAvg= Double.parseDouble(b);
//					}else if(a.equalsIgnoreCase("##mapqAvg")){
//						varMap.totalMapqAvg= Double.parseDouble(b);
//					}else if(a.equalsIgnoreCase("##readLengthAvg")){
//						varMap.readLengthAvg= Double.parseDouble(b);
//					}
//				}
//			}
//			line=bf.nextLine();
//		}
//		bf.close();
//		return varMap;
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int countNearbyVars(VarFilter varFilter) {
		return countNearbyVars(varFilter, varFilter.maxNearbyCount, varFilter.nearbyDist, 
				varFilter.nearbyGap, varFilter.flagNearby);
	}
	
	public int countNearbyVars(VarFilter varFilter, final int maxCount0, final int maxDist, final int maxGap, final boolean flag) {
		final int maxCount=maxCount0<0 ? 19 : Tools.mid(maxCount0, 8, 19);
		final Var[] array=toArray(true);
		int failed=0;
		for(int vloc=0; vloc<array.length; vloc++){
			int x=countNearbyVars(varFilter, array, vloc, maxCount, maxDist, maxGap, flag);
			if(x>maxCount){failed++;}
		}
		return failed;
	}
	
	private boolean passesSolo(Var v, VarFilter varFilter){
		assert(varFilter!=null);
		if(varFilter==null){return true;}
		boolean pass=varFilter.passesFast(v);
		if(pass){
			v.calcCoverage(scafMap);
			pass=v.forced() || varFilter.passesFilter(v, properPairRate, totalQualityAvg,
					totalMapqAvg, readLengthAvg, ploidy, scafMap, false);
		}
		return pass;
	}
	
	public int countNearbyVars(VarFilter varFilter, final Var[] array, final int vloc0, final int maxCount, 
			final int maxDist, final int maxGap, final boolean flag) {
		final Var v0=array[vloc0];
		assert(v0.nearbyVarCount==-1) : "Nearby vars were already counted?";
		int nearby=0;
		
		{//Scan left
			Var prev=v0;
			for(int i=vloc0-1; i>=0 && nearby<=maxCount; i--){
				final Var v=array[i];
				//			v.stop==v.start means adjacent;
				if(prev.start-v.stop>maxGap || v0.start-v.stop>maxDist){break;}
				
				if(!v.forced() || passesSolo(v, varFilter)){
					nearby++;
					prev=v;
				}
			}
		}
		{//Scan right
			Var prev=v0;
			for(int i=vloc0+1; i<array.length && nearby<=maxCount; i++){
				final Var v=array[i];
				//			v.stop==v.start means adjacent;
				if(v.start-prev.stop>maxGap || v.start-v0.stop>maxDist){break;}
				
				if(!v.forced() || passesSolo(v, varFilter)){
					nearby++;
					prev=v;
				}
			}
		}
		v0.nearbyVarCount=nearby;
		if(flag && nearby>varFilter.maxNearbyCount){
			v0.setFlagged(true);
		}
		return nearby;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	
	public boolean containsKey(Var v) {
		return get(v)!=null;
	}
	
	Var get(final Var v){
		final int way=v.start&MASK;
		return maps[way].get(v);
	}
	
	public long size(){
		long size=0;
		for(int i=0; i<maps.length; i++){size+=maps[i].size();}
		return size;
	}
	
	public long size2(){//123 Slow
		assert(false) : "Slow";
		int i=0;
		for(Var v : this){i++;}
		return i;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Adders            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int add(final Var v){
//		assert(size()==size2());//123;
//		assert(mappedToSelf(false));//123
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		synchronized(map){
			Var old=map.get(v);
			if(old==null){
				map.put(v, v);
//				assert(mappedToSelf(false));//123
				return 1;
			}
			else{
				synchronized(old){
					old.add(v);
//					assert(mappedToSelf(false));//123
				}
			}
		}
//		assert(size()==size2());//123;
//		assert(mappedToSelf(false));//123
		return 0;
	}
	
	int addUnsynchronized(final Var v){
//		assert(mappedToSelf(false));//123
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		Var old=map.get(v);
		if(old==null){
			map.put(v, v);
//			assert(mappedToSelf(false));//123
			return 1;
		}
		old.add(v);
//		assert(mappedToSelf(false));//123
		return 0;
	}
	
	int removeUnsynchronized(Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		return map.remove(v)==null ? 0 : 1;
	}
	
	int dumpVars(HashMap<Var, Var> mapT){
		int added=0;
		@SuppressWarnings("unchecked")
		ArrayList<Var>[] absent=new ArrayList[WAYS];
		for(int i=0; i<WAYS; i++){
			absent[i]=new ArrayList<Var>();
		}
		for(Entry<Var, Var> e : mapT.entrySet()){
			Var v=e.getValue();
			final int way=v.start&MASK;
			ConcurrentHashMap<Var, Var> map=maps[way];
			Var old=map.get(v);
			if(old==null){absent[way].add(v);}
			else{
				synchronized(old){
					old.add(v);
				}
			}
		}
		
		mapT.clear();
		for(int way=0; way<WAYS; way++){
			ConcurrentHashMap<Var, Var> map=maps[way];
			ArrayList<Var> list=absent[way];
			synchronized(map){
				for(Var v : list){
					Var old=get(v);
					if(old==null){
						map.put(v, v);
						added++;
					}
					else{
						synchronized(old){
							old.add(v);
						}
					}
				}
			}
		}
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/

	public long[] processVariantsST(VarFilter filter, long[][] scoreArray, long[] ploidyArray, long[][] avgQualityArray,
			long[] maxQualityArray, long[][] ADArray, double[] AFArray) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		long[] types=new long[Var.VAR_TYPES];
		for(ConcurrentHashMap<Var, Var> map : maps){
			long[] types2=processVariants(map, filter, null, null, null, null, null, null, false, false);
			types2=processVariants(map, filter, null, null, null, null, null, null, true, false);
			types2=processVariants(map, filter, scoreArray, ploidyArray, avgQualityArray, maxQualityArray, ADArray, AFArray, false, false);
			Tools.add(types, types2);
		}
		return types;
	}
	
//	public long[] addSharedVariantsST(VarFilter filter, VarMap sharedVarMap) {
//		assert(properPairRate>=0);
//		assert(ploidy>0);
//		assert(totalQualityAvg>=0);
//		
//		long[] types=new long[Var.VAR_TYPES];
//		for(int i=0; i<maps.length; i++){
//			long[] types2=addSharedVariants(maps[i], sharedVarMap.maps[i]);
//			Tools.add(types, types2);
//		}
//		return types;
//	}
	
	public long[] processVariantsMT(VarFilter filter, long[][] scoreArray, long[] ploidyArray, 
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray) {
		processVariantsMT_inner(filter, null, null, null, null, null, null, false);
		processVariantsMT_inner(filter, null, null, null, null, null, null, true);
		return processVariantsMT_inner(filter, scoreArray, ploidyArray, avgQualityArray, maxQualityArray, ADArray, AFArray, false);
	}
	
	private long[] processVariantsMT_inner(VarFilter filter, long[][] scoreArray, long[] ploidyArray, 
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray, boolean processInsertions) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(WAYS);
		for(int i=0; i<WAYS; i++){
			ProcessThread pt=new ProcessThread(maps[i], filter, scoreArray!=null, ploidyArray!=null, avgQualityArray!=null, ADArray!=null, processInsertions);
			alpt.add(pt);
			pt.start();
		}
		
		long[] types=new long[Var.VAR_TYPES];
		boolean success=true;
		for(ProcessThread pt : alpt){
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			if(pt.types!=null){
				Tools.add(types, pt.types);
			}
			if(scoreArray!=null){Tools.add(scoreArray, pt.scoreArray);}
			if(ploidyArray!=null){Tools.add(ploidyArray, pt.ploidyArray);}
			if(avgQualityArray!=null){Tools.add(avgQualityArray, pt.avgQualityArray);}
			if(maxQualityArray!=null){Tools.add(maxQualityArray, pt.maxQualityArray);}
			if(ADArray!=null){Tools.add(ADArray, pt.ADArray);}
			if(ADArray!=null){Tools.add(AFArray, pt.AFArray);}//Note this is triggered on ADArray
			success&=pt.success;
		}
		
		//Track whether any threads failed
//		if(!success){errorState=true;}
		
		return types;
	}
	
	private class ProcessThread extends Thread {
		
		ProcessThread(Map<Var, Var> map_, VarFilter filter_, boolean trackScores, boolean trackPloidy, 
				boolean trackQuality, boolean trackAD, boolean processInsertions_){
			map=map_;
			filter=filter_;
			scoreArray=(trackScores ? new long[8][200] : null);
			ploidyArray=(trackPloidy ? new long[ploidy+1] : null);
			avgQualityArray=(trackQuality ? new long[8][100] : null);
			maxQualityArray=(trackQuality ? new long[100] : null);
			ADArray=(trackAD ? new long[2][7] : null);
			AFArray=(trackAD ? new double[7] : null);
			processInsertions=processInsertions_;
		}
		
		@Override
		public void run(){
			types=processVariants(map, filter, scoreArray, ploidyArray, avgQualityArray, maxQualityArray, ADArray, AFArray, processInsertions, false);
			success=true;
		}
		
		final VarFilter filter;
		final Map<Var, Var> map;
		long[] types;
		final long[][] scoreArray;
		final long[] ploidyArray;
		final long[][] avgQualityArray;
		final long[] maxQualityArray;
		final long[][] ADArray;
		final double[] AFArray;
		boolean processInsertions;
		boolean success=false;
	}

	private long[] processVariants(Map<Var, Var> map, VarFilter filter, long[][] scoreArray, long[] ploidyArray, 
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray, boolean processInsertions, boolean considerNearby) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		Iterator<Entry<Var, Var>> iterator=map.entrySet().iterator();
		long[] types=new long[Var.VAR_TYPES];
		while(iterator.hasNext()){
			Entry<Var, Var> entry=iterator.next();
			final Var v=entry.getValue();
			
			if(processInsertions){
				assert(readLengthAvg>0);
				if(v.type()==Var.INS){
					synchronized(v){
						v.reviseAlleleFraction(readLengthAvg, scafMap.getScaffold(v.scafnum), this);
					}
				}
			}else{
				boolean pass=filter.passesFast(v);
				if(pass){
					v.calcCoverage(scafMap);
					pass=v.forced() || filter.passesFilter(v, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, considerNearby);
				}
				if(pass){
					types[v.type()]++;
					if(scoreArray!=null){
						int score=(int)v.phredScore(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap);
						scoreArray[0][score]++;
						scoreArray[v.type()+1][score]++;
					}
					if(ploidyArray!=null){ploidyArray[v.calcCopies(ploidy)]++;}
					if(avgQualityArray!=null){
						int q=(int)v.baseQAvg();
						avgQualityArray[0][q]++;
						avgQualityArray[v.type()+1][q]++;
					}
					if(maxQualityArray!=null){maxQualityArray[(int)v.baseQMax]++;}
					if(ADArray!=null){
						ADArray[0][v.type()]+=v.alleleCount();
						ADArray[1][v.type()]+=v.coverage();
					}
					if(AFArray!=null){AFArray[v.type()]+=v.alleleFraction();}
				}else{
					iterator.remove();
				}
			}
		}
		return types;
	}

	private long[] addSharedVariants(Map<Var, Var> map, Map<Var, Var> sharedMap) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		for(Var v : sharedMap.keySet()){
			if(!map.containsKey(v)){
				Var v2=new Var(v);
				map.put(v2, v2);
			}
		}
		
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : sharedMap.keySet()){
			v.calcCoverage(scafMap);
			types[v.type()]++;
		}
		return types;
	}
	
	public Var[] toArray(boolean sort) {
		Var[] array=new Var[(int)size()];
		int i=0;
		
//		assert(mappedToSelf(true));//123
		
		for(Var v : this){
//			System.err.println(i+"\t"+v.start+"\t"+v.stop+"\t"+v.toKey()+"\t"+new String(v.allele)+"\t"+((Object)v).hashCode());//123
			assert(i<array.length);
			array[i]=v;
			i++;
		}
		if(sort){Shared.sort(array);}
		return array;
	}
	
	private boolean mappedToSelf(boolean quiet){//123 slow
		assert(false) : "Slow";
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Var key : map.keySet()){
				Var value=map.get(key);
				assert(value!=null);
				assert(value.equals(key));
				assert(value==key);
				assert(map.get(value).equals(key));
			}
			for(Entry<Var, Var> e : map.entrySet()){
				Var key=e.getKey();
				Var value=e.getValue();
				assert(value!=null);
				assert(value.equals(key));
				assert(value==key);
			}
			for(ConcurrentHashMap<Var, Var> map2 : maps){
				if(map2!=map){
					for(Var key : map.keySet()){
						assert(!map2.containsKey(key));
					}
				}
			}
		}
		int i=0;
		for(Var v : this){
			if(!quiet){System.err.println(i+"\t"+v.start+"\t"+v.stop+"\t"+v.toKey()+"\t"+v.hashcode+"\t"+v.hashCode()+"\t"+new String(v.allele)+"\t"+((Object)v).hashCode());}
			Var v2=get(v);
			assert(v==v2);
			assert(get(v2)==v);
			assert(get(v)==v) : "\n"+i+"\t"+v2.start+"\t"+v2.stop+"\t"+v2.toKey()+"\t"+v2.hashcode+"\t"+v2.hashCode()+"\t"+new String(v2.allele)+"\t"+((Object)v2).hashCode();
			i++;
		}
		assert(i==size()) : i+", "+size()+", "+size2();
//		assert(false) : i+", "+size();
		return true;
	}
	
	public long[] calcCoverage(ScafMap scafMap) {
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : this){
			v.calcCoverage(scafMap);
			types[v.type()]++;
		}
		return types;
	}
	
	public long[] countTypes() {
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : this){
			types[v.type()]++;
		}
		return types;
	}
	
//	public void writeVarFile_(FileFormat ff, VarFilter filter, long reads, long pairs, long properPairs, long bases, String ref){
//		Var[] array=toArray(true);
//		ByteStreamWriter bsw=new ByteStreamWriter(ff);
//		bsw.start();
//		ByteBuilder bb=new ByteBuilder(33000);
//		bb.append(Var.toVarHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
//				ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
//		for(Var v : array){
//			v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap);//TODO: Track depth
//			bb.nl();
//			if(bb.length()>16384){
//				bsw.print(bb);
//				bb.clear();
//			}
//		}
//		if(bb.length()>0){
//			bsw.print(bb);
//			bb.clear();
//		}
//		bsw.poisonAndWait();
//	}
//	
//	public void writeVcfFile_(String fname, VarFilter filter, long reads, long pairs, long properPairs, long bases, String ref, String sampleName, boolean trimWhitespace){
//		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, true, false, false);
//		writeVcfFile(ff, filter, reads, pairs, properPairs, bases, ref, sampleName, trimWhitespace);
//	}
//	
//	public void writeVcfFile_(FileFormat ff, VarFilter filter, long reads, long pairs, long properPairs, long bases, String ref, String sampleName, boolean trimWhitespace){
//		Var[] array=toArray(true);
//		ByteStreamWriter bsw=new ByteStreamWriter(ff);
//		bsw.start();
//		ByteBuilder bb=new ByteBuilder(33000);
//		bb.append(Var.toVcfHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
//				ploidy, reads, pairs, properPairs, bases, ref, scafMap, sampleName, trimWhitespace)).append('\n');
//		for(Var v : array){
//			v.toVCF(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, filter, trimWhitespace);
//			bb.nl();
//			if(bb.length()>16384){
//				bsw.print(bb);
//				bb.clear();
//			}
//		}
//		if(bb.length()>0){
//			bsw.print(bb);
//			bb.clear();
//		}
//		bsw.poisonAndWait();
//	}
	
	
	public void clear() {
//		assert(mappedToSelf(false));//123
		properPairRate=-1;
		pairedInSequencingRate=-1;
		totalQualityAvg=-1;
		totalMapqAvg=-1;
		readLengthAvg=-1;
		for(int i=0; i<maps.length; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
//		assert(mappedToSelf(false));//123
	}
	
	@Override
	public String toString(){
		ByteBuilder sb=new ByteBuilder();
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Var v : map.keySet()){
				v.toTextQuick(sb);
				sb.nl();
			}
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Iteration           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public VarMapIterator iterator(){
//		if(maps.length==1){return (maps[1].entrySet().iterator();}
		return new VarMapIterator();
	}
	
	private class VarMapIterator implements Iterator<Var> {

		VarMapIterator(){
//			System.err.println("\nInit: ("+this.hashCode()+")");//123
			makeReady();
		}
		
		@Override
		public boolean hasNext() {
			return iter.hasNext();
		}

		@Override
		public Var next() {
//			System.err.println("\nNext: ("+this.hashCode()+")");//123
			Entry<Var, Var> e=iter.next();
			if(!iter.hasNext()){makeReady();}
			Var v=e==null ? null : e.getValue();
//			System.err.println("Var "+v.hashCode());//123
			return v;
		}
		
		private void makeReady(){
//			System.err.println("\nmakeReady ("+this.hashCode()+")");//123
//			System.err.println("iter="+iter+", nextMap="+nextMap);//123
			while((iter==null || !iter.hasNext()) && nextMap<maps.length){
				iter=maps[nextMap].entrySet().iterator();
				nextMap++;
//				System.err.println("iter="+iter+", nextMap="+nextMap+" (loop)");//123
			}
//			System.err.println("break");//123
		}
		
		private int nextMap=0;
		private Iterator<Entry<Var, Var>> iter=null;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public int ploidy=-1;
	public double properPairRate=-1;
	public double pairedInSequencingRate=-1;
	public double totalQualityAvg=-1;
	public double totalMapqAvg=-1;
	public double readLengthAvg=-1;
	public final ScafMap scafMap;
	final ConcurrentHashMap<Var, Var>[] maps; //ConcurrentHashMap appears to be faster than HashMap here, if there are lots of threads.
	
	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/** Must be a power of 2.  Max Vars stored is ways times 2 billion  */
	private static final int WAYS=8;
	public static final int MASK=WAYS-1;
	
}
