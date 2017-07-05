package org.ucsd.ccbb;

/**
 * SplitWorker splits a BWA file to small BAM files in the specific region.
 * 
 * @author Guorong Xu
 */
public class SplitWorker implements Runnable {
	/**a string stores the log path*/
	public String m_logFile= "";
	/**a string stores the BAM file?*/
	public String m_bamFile = "";
	/**a string stores the name of a chromosome*/
	public String m_chrom= "";
	/**an integer stores the start point of a subchromosome*/
	public int m_start = 0;
	/**an integer stores the end point of a subchromosome*/
	public int m_end = 0;
	/**a ConfigLoader stores a configuration loader*/
	public ConfigLoader m_loader = null;
	
	/**
	 * Construct a SplitWorker instance with specified configuration loader,
	 * the name of a BAM file, the name of a chromosome, the start point of a 
	 * subchromosome, the end point of a subchromosome and the log file path.
	 * @param loader	load a configuration file
	 * @param bamFile	the name of a BAM file
	 * @param chrom		a chromosome that will be splitted into equal length 
	 * 					subchromosome  
	 * @param start		the start point of a subchromosome from the chromosome
	 * @param end		the end point of a subchromosome from the chromosome
	 * @param logFile	the log file path
	 */
	public SplitWorker(ConfigLoader loader, String bamFile, String chrom, int start, int end, String logFile) {
		this.m_loader = loader;
		this.m_bamFile = bamFile;
		this.m_chrom = chrom;
		this.m_start = start;
		this.m_end = end;
		this.m_logFile = logFile;
	}
	
	@Override
	/**
	 * Overrides the function run() from interface Runnable.
	 * The function applies the S3DeleteExecutor to process commands.
	 */
    public void run() 
	{        
        //System.out.println(Thread.currentThread().getName()+ " processing file: " + m_bamFile + ".bam");
		SplitExecutor executor = new SplitExecutor();
		executor.setSAMFile(m_bamFile);
		executor.setRegion(m_chrom, m_start, m_end);
		executor.setShellFile(m_loader.get("post-alignment.split.bam.shell", "/shared/workspace/WGSPipeline/scripts/split_bam.sh"));
		executor.setLogFile(m_logFile);
		executor.execute("","");
		
		//System.out.println(Thread.currentThread().getName()+ "End.");
    }
}
