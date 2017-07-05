import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * The class is to calculate the mapping rate of alignment resulting file by novoalign for pair-end reads
 * @author guorong
 *
 */
public class CalculatePairedEndRate {
	public String m_inputFolder = "/home/workspace/data_archive/rnaSeq/VPEA/";
	public HashMap<String, ArrayList<String>> m_countsTable = new HashMap<String, ArrayList<String>>();
	
	public static void main(String[] args) throws Exception
	{
		String logFile = "/Users/Guorong/Desktop/novoalign.log";
		String suffix = ".fastq.gz";
		
		CalculatePairedEndRate parser = new CalculatePairedEndRate();
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

			System.out.println("Sample Name\tTotal Paired Reads\tPairs Aligned\tTotal Reads\tMapped Reads\t" +
					"Unique Mapped Reads\tPair Mapping Rate\tReads Mapping Rate\tUnique Mapping Rate");
			
			ArrayList<String> counts = null;
			String fastqFile = null;
			
			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("@"))
					continue;

				if (stringLine.indexOf("#  novoalignMPI -r") > -1 || stringLine.indexOf("#  novoalign -r") > -1 )
				{	
					counts = new ArrayList<String>();
					fastqFile = stringLine.substring(stringLine.indexOf(m_inputFolder) + m_inputFolder.length(), 
							stringLine.indexOf(suffix) + suffix.length());
					counts.add(fastqFile);				
				}
	
				if (stringLine.indexOf("#       Paired Reads:") > -1)
					counts.add(stringLine.substring(stringLine.indexOf("Paired Reads:") + 13));
				
				if (stringLine.indexOf("#      Pairs Aligned:") > -1)
					counts.add(stringLine.substring(stringLine.indexOf("Pairs Aligned:") + 14));
				
				if (stringLine.indexOf("#     Read Sequences:") > -1)
					counts.add(stringLine.substring(stringLine.indexOf("Read Sequences:") + 16));
				
				if (stringLine.indexOf("#            Aligned:") > -1)
					counts.add(stringLine.substring(stringLine.indexOf("Aligned:") + 9) );
				
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

