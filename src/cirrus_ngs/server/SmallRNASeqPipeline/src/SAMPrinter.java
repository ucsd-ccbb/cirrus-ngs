import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;
/**
 * The class is to print SAM file based on the sorted and normalized count tables.
 * @author guorong
 *
 */
public class SAMPrinter {
	
	public ArrayList<String[]> m_isoformKeys = new ArrayList<String[]>();
	public HashMap<String, Float> m_isomiCounts = new HashMap<String, Float>();
	public HashMap<String, String> m_samRecords = new HashMap<String, String>();
	
	public void execute(String inputPath) throws Exception
	{
		String label = "all";
		String countsType = "average";

		String countFile = inputPath + "mirRNA." + label + ".isoform_counts_" + countsType + ".txt.sort";
		String samFile = inputPath + "mirRNA." + label + ".isoforms.sam";
		String outputFile = inputPath + "mirRNA." + label + ".isoforms_counts_" + countsType + ".pileup.sam";
		
		SAMPrinter printer = new SAMPrinter();
		printer.processCountsFile(countFile.substring(0, countFile.length() - 9) + ".normalized.txt.sort");
		printer.processSAMFile(samFile);
		printer.outputSAM(outputFile);
		
		countsType = "total";
		countFile = inputPath + "mirRNA." + label + ".isoform_counts_" + countsType + ".txt.sort";
		outputFile = inputPath + "mirRNA." + label + ".isoforms_counts_" + countsType + ".pileup.sam";
		
		printer = new SAMPrinter();
		printer.processCountsFile(countFile.substring(0, countFile.length() - 9) + ".normalized.txt.sort");
		printer.processSAMFile(samFile);
		printer.outputSAM(outputFile);	
	}
	
	public void processCountsFile(String file) throws Exception
	{
		String stringLine = null;

		try
		{		
			BufferedReader in = new BufferedReader(new FileReader(file));

			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("#"))
					continue;

				StringTokenizer st = new StringTokenizer(stringLine, "\t");

				int columnth = 0;
				String miRNAName = null;
				String familyName = null;
				String isoName = null;
				String readID = null;
				String type = null;
				String counts = null;
				
				if ((stringLine.length() == 0))
					return;

				while (st.hasMoreTokens())
				{
					String stringValue = st.nextToken();

					switch (columnth)
					{
					case 0:
						miRNAName = stringValue;
						break;
					case 1:
						familyName = stringValue;
						break;		
					case 2:
						isoName = stringValue;
						break;	
					case 3:
						readID = stringValue;
						break;
					case 4:
						type = stringValue;
						break;
					case 5:
						counts = stringValue;
						break;
					default:
						;
					}
					columnth++;
				}
				
				if (familyName.indexOf("hsa-mir-548u") > -1)
					continue;
				
				m_isoformKeys.add(new String[]{miRNAName, familyName, isoName, readID, type, counts});
				
				if (!m_isomiCounts.containsKey(type + "-" + familyName))
				{
					if (Float.valueOf(counts) > 5)
						m_isomiCounts.put(type + "-" + familyName, Float.valueOf(counts));
				}
				else
				{
					if (Float.valueOf(counts) > 5)
					{
						float value = m_isomiCounts.get(type + "-" + familyName);
						m_isomiCounts.put(type + "-" + familyName, value + Float.valueOf(counts));
					}
				}
			}

			if (in != null)
				in.close();
		} catch (Exception e)
		{
			System.err.println(stringLine);
		}
	}
	
	public void processSAMFile(String file) throws Exception
	{
		String stringLine = null;

		try
		{		
			BufferedReader in = new BufferedReader(new FileReader(file));

			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("#"))
					continue;

				StringTokenizer st = new StringTokenizer(stringLine, "\t");

				int columnth = 0;
				String readID = null;
				
				if ((stringLine.length() == 0))
					return;

				while (st.hasMoreTokens())
				{
					String stringValue = st.nextToken();

					switch (columnth)
					{
					case 0:
						readID = stringValue;
						break;
					default:
						;
					}
					columnth++;
				}
				
				m_samRecords.put(readID, stringLine);
			}

			if (in != null)
				in.close();
		} catch (Exception e)
		{
			System.err.println(stringLine);
		}
	}

	public void outputSAM(String outputPath)
	{
		try
		{
			File outputFile = new File(outputPath);
			FileWriter filewriter;
			filewriter = new FileWriter(outputFile, true);
			for (int i = 0; i < m_isoformKeys.size(); i++)
			{
				//miRNAName, familyName, isoName, readID, type, counts
				//miRNAName, familyName, readID, type, counts
				
				String[] key = m_isoformKeys.get(i);
				String read = m_samRecords.get(key[3]);
				float percentage = 0;
				
				if (m_isomiCounts.containsKey(key[4] + "-" + key[1]))
					percentage =  Float.valueOf(key[5])/m_isomiCounts.get(key[4] + "-" + key[1]) * 100;
				
				DecimalFormat fnum = new DecimalFormat("##0.0");   
				
				String newRecord = read.replace(key[3], key[4] + "-(" + "num:" + key[5] + "-" + fnum.format(percentage) + "%)-" + key[0]);
				if (Float.valueOf(key[5]) > 5 && percentage > 1)
					filewriter.write( newRecord + "\n");

			}
			filewriter.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
}
