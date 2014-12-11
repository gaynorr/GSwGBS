package gbspipeline;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import net.maizegenetics.gbs.homology.TagMatchFinder;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.statistics.FisherExact;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;

/**
 *
 * @author jessepoland
 */

public class TagsToSNPsNoAnchor {
    static int minTaxaCnt=0;
    static int maxSize=100000;
    static boolean numericOutput = false;

    /**
     * @param theTBT
     * @param outHapMap
     * @param minDiff
     * @param maxMissingData
     * @param minAlleleFreq
     * @param isDHpopulation, asking if population is double haploid.  If so, will reject any SNPs that are heterozygous, 
     */
    public static void TagsToSNPsNoAnchor(TagsByTaxa theTBT, String outHapMap, int divergence, double maxMissingData, double minAlleleFreq, double maxHet, boolean isDHpopulation, boolean isBiparental, boolean callHets, double pVal) {
        
    	System.gc();
    	//do some checks
    	if(maxMissingData>1 || maxMissingData<0){
    		System.out.println("Max Missing Data (maxMissingData) must be between 0 and 1.  Entered Value: " + maxMissingData);
    		System.exit(0);
    	}
    	if(minAlleleFreq>0.5 || minAlleleFreq<0){
    		System.out.println("Minimum Allele Frequency (minAlleleFreq) must be between 0 and 0.5.  Entered Value: " + minAlleleFreq);
    		System.exit(0);
    	}
    	    	
    	final int cores=Runtime.getRuntime().availableProcessors();
        System.out.println("TagHomologyPhaseMT processors available:"+cores);
        TagMatchFinder theTMF=new TagMatchFinder(theTBT);
        ExecutorService pool = Executors.newCachedThreadPool();
        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
        FisherExact theFE=new FisherExact(theTBT.getTaxaCount()+10);
        ChiSquareTestImpl theChiSq = new ChiSquareTestImpl();

       int nSNPs=0;
       int countSNPs=0, countLoci=0;
       long time=System.currentTimeMillis();
       System.out.println("Max Pool Size "+tpe.getMaximumPoolSize());
       System.out.println("Core Pool Size "+tpe.getCorePoolSize());
       System.out.println("Comp	Tag1	Tag2	Divergence	Tag1Cnt1	Tag2Cnt2	Cnt11	exp11	Ratio11OE	P");
       
       int pauses=0;
       ArrayList<filteredSnpPair> filSnpList = new ArrayList<filteredSnpPair>(0);

       boolean[] tested = new boolean[theTBT.getTagCount()]; 
	 
        	   try {
        			BufferedWriter bw = new BufferedWriter (new FileWriter(outHapMap), 65536);
        			bw.write("rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	");
        			for (int j = 0; j < theTBT.getTaxaCount()-1; j++) {
        				bw.write(theTBT.getTaxaNames()[j]+ "\t");
        			}
        			bw.write(theTBT.getTaxaNames()[theTBT.getTaxaCount()-1]);
        			bw.newLine();
        			int count = 0;

        			int nTaxa = theTBT.getTaxaCount();
        			OpenBitSet reference = new OpenBitSet(theTBT.getTagCount());
        			
        			for (int i = 0; (i < theTBT.getTagCount())&&(nSNPs<maxSize); i++) {   
        			//	for (int i = 0; i < 100000; i++) { 
        				   long[] tagA=theTBT.getTag(i);
        				   OpenBitSet bitDistA=theTBT.getTaxaReadBitsForTag(i); // the taxa with tagA
        				   
        				   //if(bitDistA.cardinality() < nTaxa*0.3*(1-maxMissingData)) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   //if(bitDistA.cardinality() < 8) continue;  // skip the rare tags, use tags that show up at reasonable frequency (major allele)
        				   if(reference.get(i)) continue;  // check if this tag has been used as a reference, skip if so
        				   reference.set(i); // set this tags a being a reference
 
        				   TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(tagA, divergence, false);  //TreeMap is a set of Tag pairs(query Tag, hit Tag). There is only 1 mismatch between each pair
        				   if(al.size()<2) continue; // skip stuff that only has one match (self)

         				   int refCnt=0;
         				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   if(reference.get(ht.getKey())) refCnt++;
        				   }
        				   if(refCnt>1) continue;  // skip as there was already a reference tag for this TreeMap
        				   
        				   
        				   // populate an array of the tags with hits and their index
        				   ArrayList<byte[]> tags = new ArrayList<byte[]>(0);
        				   int[] hitIdx = new int[al.size()]; 
        				   //int refTagIdx = Integer.MAX_VALUE;
        				   int hi = 0;
        				   for(Entry<Integer,Integer> ht: al.entrySet()) {
        					   hitIdx[hi] = ht.getKey();
        					   hi++;
        					   long[] t = theTBT.getTag(ht.getKey());
        					   tags.add(BaseEncoder.getByteSeqFromLong(t));  
        					   //if(ht.getKey()<refTagIdx){refTagIdx = ht.getKey();} //set the first tag as the reference
        				   }
        				   //if(refTagIdx<i) continue; //skip tags that were previously tested
        				   
        				   // test each base in the set of tags for differences
        				   for(int idx=0; idx<tags.get(0).length; idx++){
        					   
        					   byte a1 = tags.get(0)[idx]; 
        					   byte a2 = tags.get(0)[idx]; // this will fail on tri-allelic sites
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]!=a1){
        							   a2 = tags.get(t)[idx];
        						   }  
        					   } 
        					   
        					   if(a1==a2) continue;  // skip over non-variable sites
        					   
        					   OpenBitSet bitDistB=theTBT.getTaxaReadBitsForTag(0);
        					   for(int z=0; z<bitDistB.size(); z++) { bitDistB.fastClear(z); }
        					   
        					   for(int t=0; t<tags.size(); t++){
        						   if(tags.get(t)[idx]==a1){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistA.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistA.fastSet(b);
        							   }
        						   }  
        						   if(tags.get(t)[idx]==a2){
        							   OpenBitSet bitDistTmp = theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int b=0; b<bitDistB.size(); b++){
        								   if(bitDistTmp.fastGet(b)) bitDistB.fastSet(b);
        							   }
        						   }
        		   
        					   } 

        		               int cntTaxa = theTBT.getTaxaCount();
        		               int cnt10 = (int)bitDistA.cardinality();
        		               int cnt01 = (int)bitDistB.cardinality();
        		               int cnt11 = (int)OpenBitSet.intersectionCount(bitDistA, bitDistB); //number of lines with both tags
        		               int cntWTag = (int)OpenBitSet.unionCount(bitDistA, bitDistB); // number of lines with one of the two tags
        		               int cnt00 = cntTaxa - cntWTag;
        					   
        					   
        		               double propPresent= (double)cntWTag/(double)cntTaxa;
        		               if(propPresent<(1-maxMissingData)) continue;  // skip if missing too much data
        		               
        		               // filter for heterozygousity
        		               if(isDHpopulation && cnt11>0) continue;  // skip if there are any heterozygous calls
        		               //if(cnt11>10) continue; // filter that the two tags don't show up in the same individuals
        		               
        		               
        		               double percentHet = (double)cnt11/(double)cntWTag;  // TODO figure out if divisor should be total taxa or #taxa with one of the tags (currently)
        		               if(percentHet>maxHet) continue; // skip if heterozygousity is > maxHet
        		            
        		               // check if the snp pair has segregation distortion, this basically does the same thing as the minor allele freq
  
        		               //if(isBiparental && (Math.abs(percent10-0.5)>0.3 || Math.abs(percent01-0.5)>0.3)) continue; // skip if one allele is >80% or <20%
        		 
        		               double freqA = (double)cnt10/(double)cntWTag; 
        		               double freqB = (double)cnt01/(double)cntWTag;
        		               if(Math.min(freqA, freqB )<minAlleleFreq) continue;
        		               
          		               double a = (double)cnt10/(double)cntTaxa;
        		               double b = (double)cnt01/(double)cntTaxa; 
        		               //double percent11 = (double)cnt11/(double)cntTaxa;
        		               //System.out.println(cnt10 + "\t" + cnt01 + "\t" + cnt11);
        		               if(cnt11 > a*b*(double)cntTaxa*0.6) continue; // skip stuff that has too much hets
        		                       		               
        		               // filter for two alleles are not associated (inbred)
        		               //double exp11=(double)cnt01*(double)cnt10/(double)theTBT.getTaxaCount();
        		               //double relExp11=(double)cnt11/exp11;
        		               double p=theFE.getCumlativeP(cnt11, cnt01, cnt10, cnt00);
        		               if(p>pVal) continue;
        		               
        						nSNPs++;
        						//String snps[] = makeSNPCalls(BaseEncoder.getSequenceFromLong(tagA), BaseEncoder.getSequenceFromLong(tagB));
        						String snps[] = new String[2];
        						snps[1] = BaseEncoder.getSequenceFromLong(tagA);             
        						snps[0] = getCharBase(a1) + "/" + getCharBase(a2);          
        						            
        						int hitIndex = hitIdx[0];
        						int div = 0;
        						
//        						filteredSnpPair fsp = new filteredSnpPair(snps[1], snps[0], i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p, bitDistA, bitDistB, callHets);
//        						fsp.swap();
//        						filSnpList.add(fsp);
        						
        						//"rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	"
        						//TODO need to check if this shoudl be  getTagLength(i)-1. That was outputting 63 bp strings
        						bw.write( BaseEncoder.getSequenceFromLong(tagA).substring(0,theTBT.getTagLength(i)) + "\t" + snps[0] + "\t" + "0" + "\t" + nSNPs + "\t" +"-1"+"\t" + idx + "\tNA\t" + cnt10 +"\t"+ cnt01 +"\t"+ cnt11 +"\t"+ propPresent +"\t" );
        						bw.write(bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOutput)+ "\n");   
        		                //System.out.println(BaseEncoder.getSequenceFromLong(tagA));
        		                
        						
        						if(nSNPs%100!=0) continue;
//        		                System.out.print(BaseEncoder.getSequenceFromLong(theTOPM.getTag(hitIndex)));
        		                System.out.println("# SNPs: " + nSNPs);
        						//System.out.printf("con %d %d %d %d %d %d %g %g %g ",i, hitIndex, div, cnt10, cnt01, cnt11, exp11, relExp11, p);
        		                //System.out.print("\n");
        		                System.out.printf("con %d %d %d %d %d %d %g ",i, hitIndex, div, cnt10, cnt01, cnt11, p);
        		                System.out.println();
        		                
        						   
        		                for(int t=0; t<tags.size(); t++){
        							   System.out.print(BaseEncoder.getSequenceFromLong(theTBT.getTag(hitIdx[t])) + "\t");
        							   OpenBitSet bitDist=theTBT.getTaxaReadBitsForTag(hitIdx[t]);
        							   for(int bt=0; bt<bitDist.size(); bt++){ System.out.print(bitDist.getBit(bt) + "\t");}
        							   System.out.println();
        						   } 
        		                System.out.println("\t\t\t\t\t\t\t\t\t" + bitsToPseudoSeq(theTBT.getTaxaCount(),bitDistA,bitDistB, snps[0], callHets, numericOutput));
        		                System.out.println(); 

        				   }
        			}

        			bw.flush();
        			bw.close();
        			System.out.println("Finished writing to: " + outHapMap);
        			System.out.println("\n# SNPs: " + nSNPs);
        			System.out.println("Divergence: " + divergence);
        			System.out.println("Het: " + maxHet);
        			System.out.println("Missing Data: " + maxMissingData);
        			System.out.println("Min Allele Freq: " + minAlleleFreq);
        			System.out.println("p-value: " + pVal);
        			System.out.println("\n" + "Finished running TagsToSNPsNoAnchor");
       
        	   }
        	   catch (Exception e) {
        			System.out.println(e.toString());
        	   }
     
		   }
    
    
    public static char getCharBase(byte base) {
        char c = 'N';
        switch (base) {
            case 0:
                return 'A';  //00
            case 1:
                return 'C';  //01
            case 2:
                return 'G';  //10
            case 3:
                return 'T';  //11
        }
        return c;
    }
    
    
	public String[] makeSNPCalls(String allele1, String allele2){
		
		String[] snpString = new String[2];
		snpString[0]="";
		snpString[1]="";
			for(int a=0; a<allele1.length(); a++){
				if(allele1.charAt(a) == allele2.charAt(a)){
					snpString[1] = snpString[1] + allele1.charAt(a);
				}
				if(!(allele1.charAt(a) == allele2.charAt(a))){
					
					snpString[0] = allele1.charAt(a) + "/" + allele2.charAt(a);
					snpString[1] = snpString[1] + "[" + allele1.charAt(a) + "/" + allele2.charAt(a) + "]";
				}
				
			}
		return snpString;
	}
    
	class filteredSnpPair implements Comparable <filteredSnpPair> {
		String sequence; //store sequence of SNP marker
		String snp;
		int queryIndex;
		int hitIndex;
		int div;
		int cCnt;
		int hCnt;
		int cnt11;
		double exp11;
		double relExp11;
		double p;
		OpenBitSet cBitDist;
		OpenBitSet hBitDist;
		
		filteredSnpPair(String sequence, String snp, int queryIndex, int hitIndex, int div, int cCnt, int hCnt, int cnt11, double exp11, double relExp11, double p, OpenBitSet cBitDist,  OpenBitSet hBitDist, boolean callHets) {
			this.sequence = sequence;
			this.snp = snp;
			this.queryIndex = queryIndex;
			this.hitIndex = hitIndex;
			this.div = div;
			this.cCnt = cCnt;
			this.hCnt = hCnt;
			this.cnt11 = cnt11;
			this.exp11 = exp11;
			this.relExp11 = relExp11;
			this.p = p;
			this.cBitDist = cBitDist;
			this.hBitDist = hBitDist;
		}
		
		public String mkStr (String name, int chr, int posi, int taxaCount, boolean callHets, boolean numericOut) {
			StringBuilder sb = new StringBuilder();
			sb.append(name).append("\t").append(snp).append("\t").append(chr).append("\t").append(posi).append("\t+\tNA\tSWGDiv\tGBS\tSWGV1\tSWGPop\tQC+\t");
			sb.append(bitsToPseudoSeq(taxaCount,cBitDist,hBitDist, snp, callHets, numericOut));
			sb.deleteCharAt(sb.length()-1);
			String str = sb.toString();
			return str;
		}
		
		public void swap () {
			if (queryIndex > hitIndex) {
				int medium;
				OpenBitSet obs;
				medium = queryIndex;
				queryIndex = hitIndex;
				hitIndex = medium;
				medium = cCnt;
				hCnt = cCnt;
				cCnt = medium;
				obs = cBitDist;
				cBitDist = hBitDist;
				hBitDist = obs;
			}
		}
		
		public int compareTo (filteredSnpPair o) {
			return queryIndex - o.queryIndex;
		}
	}
	
	public static ArrayList<filteredSnpPair> collapseSeqError(TagsByTaxa theTBT, ArrayList<filteredSnpPair> filSnpList){
		
		TagMatchFinder theTMF=new TagMatchFinder(theTBT);
		
		int i = 0;
		while(i<filSnpList.size()){
			int idx = filSnpList.get(i).queryIndex;
			long[] tagA=theTBT.getTag(idx);
			TreeMap<Integer,Integer> al=theTMF.findMatchesWithIntLengthWords(tagA, 2, false);
			OpenBitSet bitDistA=theTBT.getTaxaReadBitsForTag(i);
					
		}
		return filSnpList;
	}
	
	
    public static String bitsToPseudoSeq(int numTaxa, OpenBitSet t1, OpenBitSet t2, String snp, boolean callHets, boolean numericOut) {
        StringBuilder sb=new StringBuilder();
        //String[] allele={"A","C","R"};
        String[] tmp = snp.split("/");
        String[] allele = {tmp[0], tmp[1], "N"};
        if(callHets) allele[2] = "H";
        if(numericOut) allele = new String[]{"-1", "1", "0"};
        for (int i = 0; i < numTaxa; i++) {
            int a=t1.fastGet(i)?1:0;
            a+=t2.fastGet(i)?2:0;
            if(a==0) {
            	if(!numericOut) sb.append('N');
            	else sb.append("-9");
            		}
            //if(a==0) {sb.append("-9");}
            else {sb.append(allele[--a]);}
            sb.append("\t");
        }
        return sb.toString();
    }
    

    public static List<Object> getKeysFromValue(Map<?, ?> hm, Object value){
        List <Object>list = new ArrayList<Object>();
        for(Object o:hm.keySet()){
            if(hm.get(o).equals(value)) {
                list.add(o);
            }
        }
        return list;
    }
}
