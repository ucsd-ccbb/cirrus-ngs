
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Set;
import java.util.StringTokenizer;
/**
 * The class is to parse the gene annotation file. It creates two table, one is for all genes 
 * and the other is for all isforms. 
 * @author guorong
 *
 */
public class GFFParser2
{
	protected HashMap<String, ArrayList<String[]>> m_miRNATable = new HashMap<String, ArrayList<String[]>>();
	protected HashMap<String, String> m_miRNAIDTable = new HashMap<String, String>();
	protected HashMap<String, String[]> m_miRNATrancriptTable = new HashMap<String, String[]>();

	public void parseFile(String file) throws Exception
	{
		try
		{
			String stringLine;
			BufferedReader in = new BufferedReader(new FileReader(file));

			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("#"))
					continue;
				processGFFFile(stringLine);
			}
			
			if (in != null)
				in.close();
		} catch (Exception e)
		{
		}
	}

	private void processGFFFile(String stringLine)
	{
		StringTokenizer st = new StringTokenizer(stringLine, "\t");

		int columnth = 0;
		String chrom = null;
		String miRNAType = null;
		String start = null;
		String end = null;
		String strand = null;
		String miRNAName = null;
		String miRNAID = null;
		String primaryID = null;
		String description = null;
		
		if ((stringLine.length() == 0))
			return;

		while (st.hasMoreTokens())
		{
			String stringValue = st.nextToken();

			switch (columnth)
			{
			case 0:
				chrom = stringValue;
				break;
			case 2:
				miRNAType = stringValue;
				break;
			case 3:
				start = stringValue;
				break;
			case 4:
				end = stringValue;
				break;
			case 6:
				strand = stringValue;
				break;
			case 8:
				description = stringValue;
				break;
			default:
				;
			}
			
			columnth++;
		}
		
		if ("miRNA_primary_transcript".equalsIgnoreCase(miRNAType))
		{
			miRNAName = description.substring(description.indexOf(";Name=") + 6);
			miRNAID = description.substring(description.indexOf("ID=") + 3, description.indexOf(";accession_number"));
			
			if (!m_miRNAIDTable.containsKey(miRNAID))
			{
				String[] miRNA = new String[]{miRNAName, start, end, chrom, strand};
				m_miRNAIDTable.put(miRNAID, miRNAName);
				
				m_miRNATrancriptTable.put(miRNAName, miRNA);
				
				m_miRNATable.put(miRNAName, new ArrayList<String[]>());
			}
		}
		else 
		{
			miRNAName = description.substring(description.indexOf(";Name=") + 6, description.indexOf(";derives_from"));
			primaryID = description.substring(description.indexOf(";derives_from") + 14);
			
			if (m_miRNAIDTable.containsKey(primaryID))
			{
				String primaryMiRNAName = m_miRNAIDTable.get(primaryID);
				
				if (m_miRNATable.containsKey(primaryMiRNAName))
				{
					ArrayList<String[]> miRNAs = m_miRNATable.get(primaryMiRNAName);
					miRNAs.add(new String[]{miRNAName, start, end, chrom, strand});
				}
			}
		}
	}
	
	public String[] getPrimaryMiRNA(String miRNAName)
	{
		if (m_miRNATrancriptTable.containsKey(miRNAName))
			return m_miRNATrancriptTable.get(miRNAName);
		
		return null;
	}
	
	public String[] getMiRNA(String miRNAName, int startPosition, int endPosition) {
		String[] ret = null;
		int overlap = Integer.MIN_VALUE;
		
		if (m_miRNATable.containsKey(miRNAName))
		{
			ArrayList<String[]> miRNAs = m_miRNATable.get(miRNAName);
			
			for (int i = 0; i < miRNAs.size(); i++)
			{
				String[] miRNA = miRNAs.get(i);
				
				int start = Integer.valueOf(miRNA[1]);
				int end = Integer.valueOf(miRNA[2]);
				if ((start >= startPosition && start <= endPosition) || (end >= startPosition && end <= endPosition) 
						|| (start <= startPosition && end >= endPosition) || (start >= startPosition && end <= endPosition))
				{
					int tmpOverlap = calculateOverlap(startPosition, endPosition, start, end);
					if (tmpOverlap > overlap)
					{
						ret = miRNA;
						overlap = tmpOverlap;
					}
				}
			}
		}
		
		return ret;
	}
	
	private int calculateOverlap(int start1, int end1, int start2, int end2)
	{
		int overlap = 0;
		for (int i = start1; i < end1; i++)
			if (i >= start2 && i <= end2)
				overlap++;
		
		return overlap;
	}
	
/*	public static void main(String[] args) throws Exception {			
		String gffFile = "/Users/Guorong/Workspace/miRNA-data/hsa.gff3";
		
		GFFParser2 parser = new GFFParser2();
		parser.parseFile(gffFile);
				
		String[] mirRNA = parser.getMiRNA("hsa-mir-105-2", 151562897, 151562912);
		System.out.println(mirRNA[0] + "\t" + mirRNA[1]+ "\t" + mirRNA[2] + "\t" + mirRNA[3] + "\t" + mirRNA[4]);
	}*/
}

