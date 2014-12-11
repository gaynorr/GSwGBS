package gbspipeline;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import net.maizegenetics.gbs.pipeline.FastqToTBTPlugin;
import net.maizegenetics.gbs.pipeline.FastqToTagCountPlugin;
import net.maizegenetics.gbs.pipeline.MergeMultipleTagCountPlugin;
import net.maizegenetics.gbs.pipeline.MergeTagsByTaxaFilesPlugin;
import net.maizegenetics.gbs.pipeline.QseqToTBTPlugin;
import net.maizegenetics.gbs.pipeline.QseqToTagCountPlugin;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaUtils;

public class GBSPipeline {

    private String projectName;
    private String keyFile;
    private String seqDir;
    private String tagCountsDir;
    private String masterTagsFile;
    private String tbtDir;
    private String tbtFile;
    private String tbtFileMerge;
    private String hapDir;
    private String enzyme;
    private double maxMissingData;
    private double minorAlleleFreq;
    private double maxHet;
    private boolean isDHpopulation;
    private boolean isBiparental;
    private boolean callHets;
    private double pVal;

    public GBSPipeline(String[] rArgs) {
        projectName = rArgs[0];
        keyFile = rArgs[1];
        seqDir = rArgs[2];
        tagCountsDir = rArgs[3];
        masterTagsFile = rArgs[4];
        tbtDir = rArgs[5];
        tbtFile = rArgs[6];
        tbtFileMerge = rArgs[7];
        hapDir = rArgs[8];
        enzyme = rArgs[9];
        maxMissingData = Double.parseDouble(rArgs[10]);
        minorAlleleFreq = Double.parseDouble(rArgs[11]);
        maxHet = Double.parseDouble(rArgs[12]);
        isDHpopulation = Boolean.parseBoolean(rArgs[13]);
        isBiparental = Boolean.parseBoolean(rArgs[14]);
        callHets = Boolean.parseBoolean(rArgs[15]);
        pVal = Double.parseDouble(rArgs[16]);
    }

    public static void main(String[] args) {
        //Fastq de novo pipeline
        GBSPipeline pipeline = new GBSPipeline(args);
        pipeline.runFastqDeNovoPipeline();
    }
    
    public void runFastqDeNovoPipeline(){
        printArgs();
        runFastqToTagCountPlugin();
        runMergeMultipleTagCountPlugin();
        runFastqToTBTPlugin();
        runMergeTagsByTaxaFilesPlugin();
        mergeTaxaInTBT();
        runTagsToSNPsNoAnchor();
    }
    
    public void runQseqDeNovoPipeline(){
        printArgs();
        runQseqToTagCountPlugin();
        runMergeMultipleTagCountPlugin();
        runQseqToTBTPlugin();
        runMergeTagsByTaxaFilesPlugin();
        mergeTaxaInTBT();
        runTagsToSNPsNoAnchor(); 
    }

    public void printArgs() {
        System.out.println("\nSupplied arguments:\n");
        System.out.println("projectName " + projectName);
        System.out.println("keyFile " + keyFile);
        System.out.println("seqDir " + seqDir);
        System.out.println("tagCountsDir " + tagCountsDir);
        System.out.println("masterTagsFile " + masterTagsFile);
        System.out.println("tbtDir " + tbtDir);
        System.out.println("tbtFile " + tbtFile);
        System.out.println("tbtFileMerge " + tbtFileMerge);
        System.out.println("hapDir " + hapDir);
        System.out.println("enzyme " + enzyme);
        System.out.println("maxMissingData " + maxMissingData);
        System.out.println("minorAlleleFreq " + minorAlleleFreq);
        System.out.println("maxHet " + maxHet);
        System.out.println("isDHpopulation " + isDHpopulation);
        System.out.println("isBiparental " + isBiparental);
        System.out.println("callHets " + callHets);
        System.out.println("pVal " + pVal + "\n");
    }

    public void runFastqToTagCountPlugin() {
        String[] args = new String[]{
            "-i", seqDir,
            "-k", keyFile,
            "-e", enzyme,
            "-s", "250000000",
            "-c", "1",
            "-o", tagCountsDir,};
        FastqToTagCountPlugin plugin = new FastqToTagCountPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public void runQseqToTagCountPlugin() {
        String[] args = new String[]{
            "-i", seqDir,
            "-k", keyFile,
            "-e", enzyme,
            "-s", "250000000",
            "-c", "1",
            "-o", tagCountsDir,};
        QseqToTagCountPlugin plugin = new QseqToTagCountPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public void runMergeMultipleTagCountPlugin() {
        String[] args = new String[]{
            "-i", tagCountsDir,
            "-o", masterTagsFile,
            "-c", "16",};
        MergeMultipleTagCountPlugin plugin = new MergeMultipleTagCountPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public void runQseqToTBTPlugin() {
        String[] args = new String[]{
            "-i", seqDir,
            "-k", keyFile,
            "-e", enzyme,
            "-o", tbtDir,
            "-c", "1",
            "-t", masterTagsFile,};
        QseqToTBTPlugin plugin = new QseqToTBTPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public void runFastqToTBTPlugin() {
        String[] args = new String[]{
            "-i", seqDir,
            "-k", keyFile,
            "-e", enzyme,
            "-o", tbtDir,
            "-c", "1",
            "-t", masterTagsFile,
            "-y",};
        FastqToTBTPlugin plugin = new FastqToTBTPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public void runMergeTagsByTaxaFilesPlugin() {
        String[] args = new String[]{
            "-i", tbtDir,
            "-o", tbtFile,};
        MergeTagsByTaxaFilesPlugin plugin = new MergeTagsByTaxaFilesPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }

    public void mergeTaxaInTBT() {
        TagsByTaxaUtils.mergeTaxaByName(tbtFile,
                tbtFileMerge, FilePacking.Bit, true);
        TagsByTaxaUtils.streamBinaryToText(tbtFileMerge, 10000);
    }

    public void runTagsToSNPsNoAnchor() {
        DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd");
        Date todaysDate = new Date();
        String date = dateFormat.format(todaysDate);
        TagsByTaxa theTBT = new TagsByTaxaBitFileMap(tbtFileMerge);

        System.out.println("Starting TagsToSNPsNoAnchor nDiff=1");
        String outHapMap = hapDir + "/" + projectName + "1_" + date + ".hap";
        int nDiff = 1;
        TagsToSNPsNoAnchor.TagsToSNPsNoAnchor(theTBT, outHapMap, nDiff,
                maxMissingData, minorAlleleFreq, maxHet, isDHpopulation,
                isBiparental, callHets, pVal);
        System.gc();
        System.out.println("\n" + "DONE!");

        System.out.println("Starting TagsToSNPsNoAnchor nDiff=2");
        outHapMap = hapDir + "/" + projectName + "2_" + date + ".hap";
        nDiff = 2;
        TagsToSNPsNoAnchor.TagsToSNPsNoAnchor(theTBT, outHapMap, nDiff,
                maxMissingData, minorAlleleFreq, maxHet, isDHpopulation,
                isBiparental, callHets, pVal);
        System.gc();
        System.out.println("\n" + "DONE!");

        System.out.println("Starting TagsToSNPsNoAnchor nDiff=3");
        outHapMap = hapDir + "/" + projectName + "3_" + date + ".hap";
        nDiff = 3;
        TagsToSNPsNoAnchor.TagsToSNPsNoAnchor(theTBT, outHapMap, nDiff,
                maxMissingData, minorAlleleFreq, maxHet, isDHpopulation,
                isBiparental, callHets, pVal);
        System.gc();
        System.out.println("\n" + "DONE!");
    }
}
