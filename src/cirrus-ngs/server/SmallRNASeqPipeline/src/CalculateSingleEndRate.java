import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * The class is to calculate the mapping rate of alignment resulting file by novoalign for single-end reads
 * @author guorong
 *
 */
public class CalculateSingleEndRate {
	
	public String m_inputFolder = "/shared/workspace/software/novoIndex/hairpin_human_spike.ndx -f ";
	public HashMap<String, ArrayList<String>> m_countsTable = new HashMap<String, ArrayList<String>>();
	
	public static void main(String[] args) throws Exception
	{
		String logFile = "/Users/guorongxu/Desktop/workspace/NGSProjects/SmallRNASeq/2016_miSeqDS004_Kim_HCC_Serum_Tissue/TumorTissue_vs_NormalTissue/output/novoalign.log";
		String suffix = ".fastq";
		
		CalculateSingleEndRate parser = new CalculateSingleEndRate();
		parser.parseFile(logFile, suffix);
	}
	
	public void setInputFolder(String inputFolder)
	{
		m_inputFolder = inputFolder;
	}
	
	public void parseFile(String file, String suffix) throws Exception
	{
		try
		{
			String stringLine;
			BufferedReader in = new BufferedReader(new FileReader(file));
			
			ArrayList<String> counts = null;
			String fastqFile = null;
			
			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("@"))
					continue;
				
				if (stringLine.indexOf("#  novoalignMPI -d") > -1 || stringLine.indexOf("#  novoalign -d") > -1 )
				{	
					counts = new ArrayList<String>();
					fastqFile = stringLine.substring(stringLine.indexOf(m_inputFolder) + m_inputFolder.length(), 
							stringLine.indexOf(suffix) + suffix.length());
					counts.add(fastqFile);				
				}
				
				if (stringLine.indexOf("#     Read Sequences:") > -1)
					counts.add(stringLine.substring(stringLine.indexOf("Read Sequences:") + 16));
				
				if (stringLine.indexOf("#            Aligned:") > -1)
					counts.add(stringLine.substring(stringLine.indexOf("Aligned:") + 9));
				
				if (stringLine.indexOf("#   Unique Alignment:") > -1)
				{
					counts.add(stringLine.substring(stringLine.indexOf("Unique Alignment:") + 18));					
					m_countsTable.put(fastqFile, counts);
				}
			}

			if (in != null)
				in.close();
		} catch (Exception e)
		{
			System.err.println(e.getMessage());
		}
	}
}
