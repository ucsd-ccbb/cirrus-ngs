package org.ucsd.ccbb;

import java.util.Date;
import java.sql.Timestamp;
/**
 * WGSLauncher is the class that launches WGS pipeline with a list of BAM/FASTQ files.
 * 
 * @author	Guorong Xu
 *
 */
public class WGSLauncher {
	/**
	 * The main function calls the function process to launch WGS pipeline.
	 * The pipeline consists 8 methods of analysis
	 * @param args  the first element of the string array is the yamlFile
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception 
	{
		if (args.length < 1)
		{
			System.out.println("Usage: java -jar WGSPipelineLauncher.jar <yaml>");
			System.out.println("For example: java -jar WGSPipelineLauncher.jar /path/to/fastq2vcf.yaml");
			System.out.println("For example: java -jar WGSPipelineLauncher.jar /paty/to/bam2fastq.yaml");
			return;
		}
		String rootPath = "/shared/workspace/WGSPipeline/";
		String yamlFile = args[0];
		
		/*create a new YamlParser instance*/
		YamlParser parser = YamlParser.getInstance();
		/*parse the yaml file*/
		parser.parse(yamlFile);
		
		/*create a new WGSLauncher instance*/
		WGSLauncher launcher = new WGSLauncher();
		/*get methods of analysis*/
		String analysis = parser.getAnalysis();		
		launcher.process( rootPath,  yamlFile, analysis);
	}
	
	/**
	 * Process the yaml file by applying various methods of analysis which are listed in
	 * the yaml file
	 * @param rootPath	  the path of the configuration file
	 * @param yamlFile    a yaml file contains method of analysis, information of samples and relative project
                          information
     * @param analysis    methods of analysis in the yaml file
	 * @throws Exception  
	 */
	public void process(String rootPath, String yamlFile, String analysis) throws Exception
	{
		/*create a new ConfigLoader instance*/
		ConfigLoader loader = new ConfigLoader();
		/*process the system.conf file*/
		loader.parsePathFile(rootPath + "system.conf");
		
		/*if bam2fastq is listed in the methods of analysis,apply the 
		 *Bam2FastqCaller to convert raw BAM file to paired- end fastq files
		 */
		if ( analysis.indexOf("bam2fastq") > -1 )
		{
			Date date= new Date();
			System.out.println(new Timestamp(date.getTime()) + " System is processing Bam to Fastq...");
			Bam2FastqCaller f2bcaller = new Bam2FastqCaller();
			f2bcaller.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date.getTime()) + " Bam to Fastq is done.");
		}
		
		/*if fastqc is listed in the methods of analysis,apply the 
		 *FastQCCaller
		 */
		if ( analysis.indexOf("fastqc") > -1 )
		{
			Date date2= new Date();
			System.out.println(new Timestamp(date2.getTime()) + " System is processing fastqc...");
			FastQCCaller fcaller = new FastQCCaller();
			fcaller.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date2.getTime()) + " Fastqc is done.");
		}
		
		/*if bwa-alignment is listed in the methods of analysis,apply the 
		 *AlignmentCaller
		 */
		if ( analysis.indexOf("bwa-alignment") > -1 )
		{
			Date date3= new Date();
			System.out.println(new Timestamp(date3.getTime()) + " System is processing bwa alignment...");
			AlignmentCaller bwa = new AlignmentCaller();
			bwa.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date3.getTime()) + " Bwa alignment is done.");
		}
		
		/*if post-alignment is listed in the methods of analysis,apply the 
		 *PostAlignCaller
		 */
		if ( analysis.indexOf("post-alignment") > -1 )
		{
			Date date4= new Date();
			System.out.println(new Timestamp(date4.getTime()) + " System is processing post-alignment...");
			PostAlignCaller bcaller = new PostAlignCaller();
			bcaller.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date4.getTime()) + " Post-alignment is done.");
		}
		
		/*if gatk-haplotype is listed in the methods of analysis,apply the 
		 *HaplotyperCaller to use haplotype to call variants
		 */
		if ( analysis.indexOf("gatk-haplotype") > -1 )
		{
			Date date5= new Date();
			System.out.println(new Timestamp(date5.getTime()) + " System is processing GATK-haplotype...");
			HaplotypeCaller hcaller = new HaplotypeCaller();
			hcaller.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date5.getTime()) + " GATK-haplotype is done.");
		}
		
		/*if both of gatk-haplotype and merge is listed in the methods of 
		 *analysis,apply the MergeCaller 
		 */
		if ( analysis.indexOf("gatk-haplotype") > -1 || analysis.indexOf("merge") > -1 )
		{
			Date date6= new Date();
			System.out.println(new Timestamp(date6.getTime()) + " System is processing merge BAM and VCF files...");
			MergeCaller mcaller = new MergeCaller();
			mcaller.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date6.getTime()) + " Merging BAM and VCF files is done.");
		}
		
		/*if combine-vcf is listed in the methods of analysis,apply the 
		 *GroupSampleCaller to call gatk CombineGVCFs to group VCF files into
		 *one group.
		 */
		if ( analysis.indexOf("combine-vcf") > -1 )
		{
			Date date7= new Date();
			System.out.println(new Timestamp(date7.getTime()) + " System is grouping VCF files...");
			GroupSampleCaller gcaller = new GroupSampleCaller();
			gcaller.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date7.getTime()) + " Grouping VCF files is done.");
		}
		
		/*if variant-filter is listed in the methods of analysis,apply the 
		 *VariantFilterCaller to call gatk VariantRecalibrator to filter 
		 *variants
		 */
		if ( analysis.indexOf("variant-filter") > -1 )
		{
			Date date8= new Date();
			System.out.println(new Timestamp(date8.getTime()) + " System is filtering VCF files...");
			VariantFilterCaller vfcaller = new VariantFilterCaller();
			vfcaller.run(rootPath, loader, yamlFile);
			System.out.println(new Timestamp(date8.getTime()) + " Filtering VCF files is done.");
		}
	}
}