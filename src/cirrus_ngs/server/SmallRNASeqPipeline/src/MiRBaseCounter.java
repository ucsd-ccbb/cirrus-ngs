import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.StringTokenizer;
/**
 * The class is to process the alignment resulting files against complete mirbase reference.
 * @author guorong
 *
 */
public class MiRBaseCounter {
	public ArrayList<String> m_fileList = new ArrayList<String>();
	public ArrayList<String> m_fileNameList = new ArrayList<String>();
	public ArrayList<String> m_miRNAList = new ArrayList<String>();
	public ArrayList<String> m_mappableReadsList = new ArrayList<String>();
	public HashMap<String, ArrayList<Integer>> m_miRNACountsList = new HashMap<String, ArrayList<Integer>>();
	public HashMap<String, ArrayList<String>> m_miRNAAnnotationList = new HashMap<String, ArrayList<String>>();
	public HashMap<String, String> m_miRNAReadsList = new HashMap<String, String>();
	public HashMap<String, String> m_miRNAReadsIDList = new HashMap<String, String>();
	public boolean m_ignoreMismatchLocation = false;
	public MiRNAUtil m_miRNAUtil = MiRNAUtil.getInstance();
	
	public GFFParser2 m_gffParser = null;
	public int m_totalReads = 0;
	public int m_totalUnmappedReads = 0;
	public int m_totalIgnoredReads = 0;
	public int m_totalIgnoredLess17 = 0;
	public int m_totalIgnoredLonger27 = 0;
	public int m_totalMismatch2 = 0;
	
	public static void main(String[] args) throws Exception {		
		String gffFile = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/data/miRNASeq/hsa.gff3";
		String file = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/data/miRNASeq/AD5-WK24_GTGGCC_S67_L002_R1_001.trim.fastq.sam";
		int fileIndex = 0;
		
		GFFParser2 parser = new GFFParser2();
		parser.parseFile(gffFile);
		
		MiRBaseCounter counter = new MiRBaseCounter();
		counter.setGFFParser(parser);
		counter.m_fileNameList.add("AD5-WK24_GTGGCC_S67_L002_R1_001.trim.fastq.sam");
		counter.parseFile(file,  fileIndex);
	}
	
	public void setGFFParser(GFFParser2 parser)
	{
		m_gffParser = parser;
	}
	
	public void setIgnoreMismatchLocation(boolean ignoreMismatchLocation)
	{
		m_ignoreMismatchLocation = ignoreMismatchLocation;
	}
	
	public void listAllFiles(String folder, String suffix)
	{
		File[] listOfFiles = new File(folder).listFiles();
		if (listOfFiles == null || listOfFiles.length == 0)
			return;

		for (int i = 0; i < listOfFiles.length; i++)
		{
			if (listOfFiles[i].isFile())
			{
				String fileName = listOfFiles[i].getName();

				if (fileName.endsWith(suffix))
				{
					m_fileNameList.add(fileName);
					m_fileList.add(listOfFiles[i].getPath());
					//System.out.println(listOfFiles[i].getPath());
				}
			} else if (listOfFiles[i].isDirectory())
			{
				listAllFiles(listOfFiles[i].getPath(), suffix);
			}
		}
	}

	public void parseFile(String file, int fileIndex) throws Exception
	{
		String stringLine = null;
		m_totalReads = 0;
		m_totalUnmappedReads = 0;
		m_totalIgnoredReads = 0;
		m_totalIgnoredLess17 = 0;
		m_totalIgnoredLonger27 = 0;
		m_totalMismatch2 = 0;
		
		try
		{		
			BufferedReader in = new BufferedReader(new FileReader(file));

			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("@"))
					continue;

				processSAMFile(stringLine, fileIndex);
			}

			if (in != null)
				in.close();
		} catch (Exception e)
		{
			System.err.println(stringLine);
		}
	}

	private void processSAMFile(String stringLine, int fileIndex)
	{
		StringTokenizer st = new StringTokenizer(stringLine, "\t");

		int columnth = 0;
		String miRNAName = null;
		String flag = null;
		int position = 0;
		int readLength = 0;
		String rawSequence = null;
		int matchedLength = 0;
		String readID = null;
		String cigar = null;
		String MDString = null;
		
		if ((stringLine.length() == 0))
			return;

		while (st.hasMoreTokens())
		{
			String stringValue = st.nextToken();

			switch (columnth)
			{
			case 0:
				readID = stringValue;
				
				m_totalReads++;
				break;
			case 1:
				flag = stringValue;
				if (!"0".equalsIgnoreCase(flag) && !"16".equalsIgnoreCase(flag))
				{
					m_totalUnmappedReads++;
					return;
				}
				
				break;			
			case 2:
				miRNAName = stringValue;
				if ("*".equalsIgnoreCase(miRNAName))
				{
					m_totalUnmappedReads++;
					return;
				}
				break;
			case 3:
				position = Integer.valueOf(stringValue);
				break;
			case 5:
				cigar = stringValue;
				int location = stringValue.indexOf("M");
				String lengthString  =  stringValue.substring(location - 2, location);
				matchedLength = Integer.valueOf(lengthString);
				break;
			case 9:
				rawSequence = stringValue;
				readLength = stringValue.length();
				break;
			case 12:
				MDString = stringValue;
				break;
			default:
				;
			}
			columnth++;
		}
		
		int mismatchNum = MiRNAUtil.getInstance().getMismatchNum(MDString);

		if (mismatchNum > 1)
		{
			m_totalMismatch2++;
		}
		
		int length = readLength;
		
		if (length < 17)
		{
			m_totalIgnoredLess17++;
			return;
		}
		
		if (length > 27)
		{
			m_totalIgnoredLonger27++;
			return;
		}
		
		if (matchedLength > 0 && matchedLength < readLength)
			length = matchedLength;
		
		String[] miRNA = m_gffParser.getPrimaryMiRNA(miRNAName);
		
		if (miRNA == null)
		{
			m_totalIgnoredReads++;
			if (!m_miRNACountsList.containsKey(miRNAName))
			{
				ArrayList<Integer> countsList = new ArrayList<Integer>();
				for (int i = 0; i < m_fileNameList.size(); i++)
					countsList.add(0);
				
				countsList.set(fileIndex, 1);
				
				m_miRNACountsList.put(miRNAName, countsList);
				
				ArrayList<String> ann = new ArrayList<String>();
				ann.add(miRNAName);
				ann.add(miRNAName);
				ann.add(miRNAName);
				ann.add(miRNAName);
				ann.add(String.valueOf(position - 1));
				ann.add(String.valueOf(position + length - 2));
				ann.add(String.valueOf(0));
				ann.add(String.valueOf(0));
				ann.add(rawSequence);
				m_miRNAAnnotationList.put(miRNAName, ann);
			}
			else
			{
				ArrayList<Integer> countsList = m_miRNACountsList.get(miRNAName);
				countsList.set(fileIndex, countsList.get(fileIndex) + 1);	
			}
			
			if (!m_miRNAReadsList.containsKey(miRNAName))
			{
				m_miRNAReadsList.put(miRNAName, stringLine);
				m_miRNAReadsIDList.put(miRNAName, readID);
			}
			
			return;
		}
		
		// miRNAName, start, end, chrom, strand
		String[] miRNAIsoform = m_gffParser.getMiRNA(miRNAName, Integer.valueOf(miRNA[1]) + position - 1, 
				Integer.valueOf(miRNA[1]) + position + length - 2);
		
		if ("-".equalsIgnoreCase(miRNA[4]))
			miRNAIsoform = m_gffParser.getMiRNA(miRNAName, Integer.valueOf(miRNA[2]) - position - length - 2, 
				Integer.valueOf(miRNA[2]) - position - 1);
		
		if (miRNAIsoform != null)
		{
			String miRNAId = miRNAName + ":" + miRNAIsoform[0];
			int left = 0;
			int right = 0;
			
			if ("+".equalsIgnoreCase(miRNA[4]))
			{
				left = Integer.valueOf(miRNA[1]) + position - Integer.valueOf(miRNAIsoform[1]) - 1;
				if (left >= 0)
					miRNAId = miRNAId + "-left+" + left;
				else if (left < 0)
					miRNAId = miRNAId + "-left" + left;
				
				right = Integer.valueOf(miRNA[1]) + position + length - Integer.valueOf(miRNAIsoform[2]) - 2;
				
				if (right >= 0)
					miRNAId = miRNAId + "-right+" + right;
				else if (right < 0)
					miRNAId = miRNAId + "-right" + right;
			}
			else
			{
				left = Integer.valueOf(miRNA[2]) - position - Integer.valueOf(miRNAIsoform[2]) + 1;
				if (left >= 0)
					miRNAId = miRNAId + "-left-" + left;
				else if (left < 0)
					miRNAId = miRNAId + "-left+" + Math.abs(left);
				
				right = Integer.valueOf(miRNA[2]) - position - length - Integer.valueOf(miRNAIsoform[1]) + 2;
				if (right >= 0)
					miRNAId = miRNAId + "-right-" + right;
				else if (right < 0)
					miRNAId = miRNAId + "-right+" + Math.abs(right);
			}
			
			if (left == 0 && right == 0)
				miRNAId = miRNAName + ":" + miRNAIsoform[0];
			
			//Find out the mismatch location and create a new isoform record.
			
			if (mismatchNum > 0)
			{
				Integer[] locations = MiRNAUtil.getInstance().getAbsoluteMismatchLocation(MDString);
				for (int k = 0; k < locations.length; k++)
				{
					int absoluteLocation = position + locations[k];					
					int relativeLocation = MiRNAUtil.getInstance().getRelativeMirBaseMismatchLocation(absoluteLocation, miRNA, miRNAIsoform);
					
					if (!m_ignoreMismatchLocation && Math.abs(relativeLocation) > 8)
						return;
					
					if (relativeLocation > 0)
						miRNAId = miRNAId + ":M+" + relativeLocation;
					else
						miRNAId = miRNAId + ":M" + (relativeLocation -1);
				}
			}
			
			if (!m_miRNACountsList.containsKey(miRNAId))
			{
				ArrayList<Integer> countsList = new ArrayList<Integer>();
				for (int i = 0; i < m_fileNameList.size(); i++)
					countsList.add(0);
				
				countsList.set(fileIndex, 1);
				
				m_miRNACountsList.put(miRNAId, countsList);
				
				ArrayList<String> ann = new ArrayList<String>();
				ann.add(miRNA[0]);
				ann.add(miRNAIsoform[0]);
				ann.add(miRNA[3]);
				ann.add(miRNA[4]);
				ann.add(miRNAIsoform[1]);
				ann.add(miRNAIsoform[2]);
				if ("+".equalsIgnoreCase(miRNA[4]))
				{
					ann.add(String.valueOf(left));
					ann.add(String.valueOf(right));
				}
				else
				{
					ann.add(String.valueOf(-left));
					ann.add(String.valueOf(-right));
				}
				ann.add(rawSequence);
				
				m_miRNAAnnotationList.put(miRNAId, ann);
				if (!m_miRNAReadsList.containsKey(miRNAId))
				{
					m_miRNAReadsList.put(miRNAId, stringLine);
					m_miRNAReadsIDList.put(miRNAId, readID);
				}
			}
			else 
			{
				ArrayList<Integer> countsList = m_miRNACountsList.get(miRNAId);
				countsList.set(fileIndex, countsList.get(fileIndex) + 1);
			}
		}
		else
		{
			m_totalIgnoredReads++;			
			
			if (!m_miRNACountsList.containsKey(miRNAName + ":undefined"))
			{
				ArrayList<Integer> countsList = new ArrayList<Integer>();
				for (int i = 0; i < m_fileNameList.size(); i++)
					countsList.add(0);
				
				countsList.set(fileIndex, 1);
				
				m_miRNACountsList.put(miRNAName + ":undefined", countsList);
				
				ArrayList<String> ann = new ArrayList<String>();
				ann.add(miRNA[0]);
				ann.add(miRNA[0]);
				ann.add(miRNA[3]);
				ann.add(miRNA[4]);
				if ("+".equalsIgnoreCase(miRNA[4]))
				{
					ann.add(String.valueOf(Integer.valueOf(miRNA[1]) + position - 1));
					ann.add(String.valueOf(Integer.valueOf(miRNA[1]) + position + length - 2));
				}
				else
				{
					ann.add(String.valueOf(Integer.valueOf(miRNA[2]) - position + 1));
					ann.add(String.valueOf(Integer.valueOf(miRNA[2]) - position - length + 2));
				}
				ann.add(String.valueOf(0));
				ann.add(String.valueOf(0));
				ann.add(rawSequence);
				m_miRNAAnnotationList.put(miRNAName + ":undefined", ann);
				
				if (!m_miRNAReadsList.containsKey(miRNAName + ":undefined"))
				{
					m_miRNAReadsList.put(miRNAName + ":undefined", stringLine);
					m_miRNAReadsIDList.put(miRNAName + ":undefined", readID);
				}
			}
			else
			{
				ArrayList<Integer> countsList = m_miRNACountsList.get(miRNAName + ":undefined");
				countsList.set(fileIndex, countsList.get(fileIndex) + 1);	
			}
			return;		
		}
	}
	
	public void loadDesignTable(String designFile)
	{	
		try {
			YamlParser parser = YamlParser.getInstance();
			m_miRNAUtil.setGroupList(parser.m_sampleGroups);
			m_miRNAUtil.setGroupTypes(parser.m_groupTypes);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
