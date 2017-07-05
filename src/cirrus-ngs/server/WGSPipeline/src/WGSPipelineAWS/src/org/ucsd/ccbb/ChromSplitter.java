package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * ChromSplitter splits a chromosome into subchromosome.
 * 
 * @author Guorong Xu
 */
public class ChromSplitter {
	/**an arraylist of string arrays stores equal length subchromosome on each
	 * chromosome */
	public ArrayList<String[]> m_regions = null;
	/**an arraylist of strings stores names of 25 chromosomes*/
	public ArrayList<String> m_chroms = null;
	/**a array of string arrays stores names and lengths of 25 chromosomes*/
	public String[][] m_chromLength = new String[][]{
			{"chrM", "16571"},
			{"chr1", "249250621"},
			{"chr2", "243199373"},
			{"chr3", "198022430"},
			{"chr4", "191154276"},
			{"chr5", "180915260"},
			{"chr6", "171115067"},
			{"chr7", "159138663"},
			{"chr8", "146364022"},
			{"chr9", "141213431"},
			{"chr10", "135534747"},
			{"chr11", "135006516"},
			{"chr12", "133851895"},
			{"chr13", "115169878"},
			{"chr14", "107349540"},
			{"chr15", "102531392"},
			{"chr16", "90354753"},
			{"chr17", "81195210"},
			{"chr18", "78077248"},
			{"chr19", "59128983"},
			{"chr20", "63025520"},
			{"chr21", "48129895"},
			{"chr22", "51304566"},
			{"chrX", "155270560"},
			{"chrY", "59373566"},
	};

	/**
	 * This function split each chromosome into specified equal length 
	 * subchromosomes.
	 * @param loader	load the configuration file
	 * @return			an arraylist of string arrays storing equal length 
	 * 					subchromosome on each chromosome 
	 * @throws NumberFormatException	an exception will be thrown if trying 
	 * 									to convert a string to integer, but the
	 * 									string does not have the appropriate 
	 * 									format
	 * @throws InterruptedException		an exception will be thrown when a 
	 * 									thread is waiting, sleeping, or 
	 * 									otherwise occupied, and the thread is 
	 * 									interrupted, either before or during 
	 * 									the activity
	 */
	public ArrayList<String[]> getRegions(ConfigLoader loader) throws NumberFormatException, InterruptedException 
	{
		m_regions = new ArrayList<String[]>();
		int interval = Integer.valueOf(loader.get("sam.chrom.length.interval", "30000000"));
		
		/*m_chromLength.length is 25*/
		for (int i = 0; i < m_chromLength.length; i ++)
		{
			/*get the name of the chromosome*/
			String chrom = m_chromLength[i][0];
			/*get the length of the chromosome*/
			int length = Integer.valueOf(m_chromLength[i][1]);
			int regionNum = length / interval + 1;
			
			for (int j = 0; j < regionNum; j ++)
			{
				int start = j * interval + 1;
				int end = (j + 1) * interval;
				
				if (end > length)
					end = length;
				
				m_regions.add(new String[]{chrom, String.valueOf(start), String.valueOf(end)});
			}
		}
		return m_regions;
	}
	
	/**
	 * This function gets names of 25 chromosomes.
	 * @return	an arraylist of strings storing names of 25 chromosomes
	 * @throws NumberFormatException
	 * @throws InterruptedException
	 */
	public ArrayList<String> getChoms() throws NumberFormatException, InterruptedException 
	{
		m_chroms = new ArrayList<String>();
		
		for (int i = 0; i < m_chromLength.length; i ++)
		{
			m_chroms.add(m_chromLength[i][0]);
		}
		
		return m_chroms;
	}
	
}
