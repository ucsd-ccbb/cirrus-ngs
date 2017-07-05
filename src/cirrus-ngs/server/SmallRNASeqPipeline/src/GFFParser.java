import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Set;
import java.util.StringTokenizer;
/**
 * The class is to parse the gene annotation file. It creates two table, one is for all genes on the forward strand 
 * and the other is for all genes on the reverse strand. 
 * @author guorong
 *
 */
public class GFFParser
{
	protected HashMap<String, ArrayList<String[]>> m_miRNAForwardTable = new HashMap<String, ArrayList<String[]>>();
	protected HashMap<String, ArrayList<String[]>> m_miRNAReverseTable = new HashMap<String, ArrayList<String[]>>();
	protected HashMap<String, String[]> m_miRNATranHash = new HashMap<String, String[]>();
	protected HashMap<String, ArrayList<String[]>> m_miRNAForwardTranList = new HashMap<String, ArrayList<String[]>>();
	protected HashMap<String, ArrayList<String[]>> m_miRNAReverseTranList = new HashMap<String, ArrayList<String[]>>();
	
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
			
			sort(m_miRNAForwardTable);
			sort(m_miRNAReverseTable);
			
			sort(m_miRNAForwardTranList);
			sort(m_miRNAReverseTranList);
			
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
				miRNAName = description.substring(description.indexOf(";Name=") + 6);
				if (miRNAName.indexOf(";derives_from") > -1)
					miRNAName = miRNAName.substring(0, miRNAName.indexOf(";derives_from"));
				break;
			default:
				;
			}
			
			columnth++;
		}
		
		if ("miRNA_primary_transcript".equalsIgnoreCase(miRNAType))
		{
			if (!m_miRNATranHash.containsKey(miRNAName))
			{
				String[] miRNA = new String[]{miRNAName, start, end, chrom, strand};
				m_miRNATranHash.put(miRNAName, miRNA);
				
				if ("+".equalsIgnoreCase(strand))
				{
					if (m_miRNAForwardTranList.containsKey(chrom))
					{
						ArrayList<String[]> miRNAList = m_miRNAForwardTranList.get(chrom);
						String[] miRNATran = new String[]{miRNAName, start, end, chrom, strand};
						miRNAList.add(miRNATran);
					}
					else 
					{
						ArrayList<String[]> miRNAList = new ArrayList<String[]>();
						String[] miRNATran = new String[]{miRNAName, start, end, chrom, strand};
						miRNAList.add(miRNATran);
						
						m_miRNAForwardTranList.put(chrom, miRNAList);
					}
				}
				else if ("-".equalsIgnoreCase(strand))
				{
					if (m_miRNAReverseTranList.containsKey(chrom))
					{
						ArrayList<String[]> miRNAList = m_miRNAReverseTranList.get(chrom);
						String[] miRNATran = new String[]{miRNAName, start, end, chrom, strand};
						miRNAList.add(miRNATran);
					}
					else 
					{
						ArrayList<String[]> miRNAList = new ArrayList<String[]>();
						String[] miRNATran = new String[]{miRNAName, start, end, chrom, strand};
						miRNAList.add(miRNATran);
						
						m_miRNAReverseTranList.put(chrom, miRNAList);
					}
				}
			}
		}
		else 
		{
			if ("+".equalsIgnoreCase(strand))
			{
				if (m_miRNAForwardTable.containsKey(chrom))
				{
					ArrayList<String[]> miRNAList = m_miRNAForwardTable.get(chrom);
					String[] miRNA = new String[]{miRNAName, start, end};
					miRNAList.add(miRNA);
				}
				else 
				{
					ArrayList<String[]> miRNAList = new ArrayList<String[]>();
					String[] miRNA = new String[]{miRNAName, start, end};
					miRNAList.add(miRNA);
					
					m_miRNAForwardTable.put(chrom, miRNAList);
				}
			}
			else if ("-".equalsIgnoreCase(strand))
			{
				if (m_miRNAReverseTable.containsKey(chrom))
				{
					ArrayList<String[]> miRNAList = m_miRNAReverseTable.get(chrom);
					String[] miRNA = new String[]{miRNAName, start, end, chrom, strand};
					miRNAList.add(miRNA);
				}
				else 
				{
					ArrayList<String[]> miRNAList = new ArrayList<String[]>();
					String[] miRNA = new String[]{miRNAName, start, end, chrom, strand};
					miRNAList.add(miRNA);
					
					m_miRNAReverseTable.put(chrom, miRNAList);
				}
			}
		}
	}
	
	public String[] getMiRNA(String miRNAName)
	{
		if (m_miRNATranHash.containsKey(miRNAName))
				return m_miRNATranHash.get(miRNAName);
		
		return null;
	}
	
	public String[] getPrimaryMiRNA(String chrom, String strand, String start, String end)
	{
		String[] miRNA = null;
		if ("+".equalsIgnoreCase(strand) && m_miRNAForwardTranList.containsKey(chrom))
			miRNA = lookforMiRNA(m_miRNAForwardTranList.get(chrom), Integer.valueOf(start), Integer.valueOf(end));
		else if ("-".equalsIgnoreCase(strand) && m_miRNAReverseTranList.containsKey(chrom))
			miRNA = lookforMiRNA(m_miRNAReverseTranList.get(chrom), Integer.valueOf(start), Integer.valueOf(end));

		return miRNA;
	}
	
	public String[] getMiRNA(String chrom, String strand, String start, String end)
	{
		String[] miRNA = null;
		if ("+".equalsIgnoreCase(strand) && m_miRNAForwardTable.containsKey(chrom))
			miRNA = lookforMiRNA(m_miRNAForwardTable.get(chrom), Integer.valueOf(start), Integer.valueOf(end));
		else if ("-".equalsIgnoreCase(strand) && m_miRNAReverseTable.containsKey(chrom))
			miRNA = lookforMiRNA(m_miRNAReverseTable.get(chrom), Integer.valueOf(start), Integer.valueOf(end));

		return miRNA;
	}
	
	private String[] lookforMiRNA(ArrayList<String[]> miRNAs, int startPosition, int endPosition) {
		int low = 0;
		int high = miRNAs.size() - 1;
		int start = Integer.valueOf(miRNAs.get(low)[1]);
		int end = Integer.valueOf(miRNAs.get(low)[2]);
		if ((start >= startPosition && start <= endPosition) || (end >= startPosition && end <= endPosition) 
				|| (start <= startPosition && end >= endPosition) || (start >= startPosition && end <= endPosition))
			return miRNAs.get(low);
		
		start = Integer.valueOf(miRNAs.get(high)[1]);
		end = Integer.valueOf(miRNAs.get(high)[2]);
		if ((start >= startPosition && start <= endPosition) || (end >= startPosition && end <= endPosition) 
				|| (start <= startPosition && end >= endPosition) || (start >= startPosition && end <= endPosition))
			return miRNAs.get(high);
		
		while (high >= low) {
			int middle = (low + high) / 2;
			start = Integer.valueOf(miRNAs.get(middle)[1]);
			end = Integer.valueOf(miRNAs.get(middle)[2]);

			if (end < startPosition) 
				low = middle + 1;
			else if (start > endPosition) 
				high = middle - 1;
			else if ((start >= startPosition && start <= endPosition) || (end >= startPosition && end <= endPosition) 
					|| (start <= startPosition && end >= endPosition) || (start >= startPosition && end <= endPosition))
				return miRNAs.get(middle);
		}
		return null;
	}
	
	private void sort(HashMap<String, ArrayList<String[]>> mirRNATable) {
		Set<String> keys = mirRNATable.keySet();
		// Iterate over the keys
		for (String key : keys) {
			ArrayList<String[]> mirRNAs = mirRNATable.get(key);
			Collections.sort(mirRNAs, new Comparator<String[]>() {
				public int compare(String[] strings, String[] otherStrings) {
					return Integer.valueOf(strings[1]).compareTo(Integer.valueOf(otherStrings[1]));
				}
			});
		}
	}
	
	public static void main(String[] args) throws Exception {			
		String gffFile = "/Users/Guorong/Workspace/miRNA-data/hsa.gff3";
		
		GFFParser parser = new GFFParser();
		parser.parseFile(gffFile);
				
		String[] mirRNA = parser.getPrimaryMiRNA("chr10", "-", "97824126", "97824147");
		System.out.println(mirRNA[0] + "\t" + mirRNA[1]+ "\t" + mirRNA[2]);
		
		mirRNA = parser.getMiRNA("chr10", "-", "97824126", "97824147");
		System.out.println(mirRNA[0] + "\t" + mirRNA[1]+ "\t" + mirRNA[2]);
	}
}

