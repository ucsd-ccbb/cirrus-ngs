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
 * The class is to process the alignment resulting files against complete genome reference.
 * @author guorong
 *
 */
public class GenomeCounter {
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
	
	public GFFParser m_gffParser = null;
	public int m_totalReads = 0;
	public int m_totalUnmappedReads = 0;
	public int m_totalIgnoredReads = 0;
	public int m_totalIgnoredLess17 = 0;
	public int m_totalIgnoredLonger27 = 0;
	public int m_totalMismatch2 = 0;
	
	public void setGFFParser(GFFParser parser)
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

				if (fileName.endsWith(suffix) && m_miRNAUtil.m_groupList.containsKey(fileName))
				{
					m_fileNameList.add(fileName);
					m_fileList.add(listOfFiles[i].getPath());
					System.out.println(listOfFiles[i].getPath());
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

	// The method is to process the alignment resulting files by novoalign. 
	// And it also calculate the coordinator shift according to the gene annotation file.
	private void processSAMFile(String stringLine, int fileIndex)
	{
		StringTokenizer st = new StringTokenizer(stringLine, "\t");

		int columnth = 0;
		String chrom = null;
		String flag = null;
		int position = 0;
		int readLength = 0;
		int matchedLength = 0;
		String rawSequence = null;
		String readID = null;
		String cigar = null;
		String NMString = null;
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
				chrom = stringValue;
				if ("*".equalsIgnoreCase(chrom))
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
				int lastlocation = stringValue.lastIndexOf("M");
				if (lastlocation > -1 && Character.isDigit(lastlocation - 2))
				{
					String lengthString  =  stringValue.substring(lastlocation - 2, lastlocation);
					matchedLength = Integer.valueOf(lengthString);
				}
				
				if (location > -1 && Character.isDigit(lastlocation - 2))
				{
					String lengthString  =  stringValue.substring(location - 2, location);
					matchedLength = Integer.valueOf(lengthString);					
				}
				break;
			case 9:
				rawSequence = stringValue;
				readLength = stringValue.length();
				break;
			case 14:
				NMString = stringValue;
				break;
			case 15:
				MDString = stringValue;
				break;
			default:
				;
			}
			columnth++;
		}
		
		int mismatchNum = MiRNAUtil.getInstance().getMismatchNum(NMString);
		if (mismatchNum > 1)
			m_totalMismatch2++;
		
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
		
		String[] miRNA = null;
		if ("0".equalsIgnoreCase(flag) )
			miRNA = m_gffParser.getPrimaryMiRNA(chrom, "+", String.valueOf(position), String.valueOf(position + length));
		else if ("16".equalsIgnoreCase(flag))
			miRNA = m_gffParser.getPrimaryMiRNA(chrom, "-", String.valueOf(position), String.valueOf(position + length));

		if (miRNA == null)
		{
			m_totalIgnoredReads++;
			return;
		}
		
		// miRNAName, start, end, chrom, strand
		String[] miRNAIsoform = m_gffParser.getMiRNA(chrom, "+", String.valueOf(position), 
				String.valueOf(position + length - 1));	

		if ("-".equalsIgnoreCase(miRNA[4]))
			miRNAIsoform = m_gffParser.getMiRNA(chrom, "-",String.valueOf(position), 
					String.valueOf(position + length - 1));
			
		if (miRNAIsoform != null)
		{
			String miRNAId = miRNA[0] + ":" + miRNAIsoform[0];
			int left = 0;
			int right = 0;
			
			if ("+".equalsIgnoreCase(miRNA[4]))
			{
				left = position - Integer.valueOf(miRNAIsoform[1]);
				if (left >= 0)
					miRNAId = miRNAId + "-left+" + left;
				else if (left < 0)
					miRNAId = miRNAId + "-left" + left;
				
				right = position + length - Integer.valueOf(miRNAIsoform[2]) - 1;
				
				if (right >= 0)
					miRNAId = miRNAId + "-right+" + right;
				else if (right < 0)
					miRNAId = miRNAId + "-right" + right;
			}
			else
			{
				left = position + length - Integer.valueOf(miRNAIsoform[2]) - 1;
				
				if (left >= 0)
					miRNAId = miRNAId + "-left-" + left;
				else if (left < 0)
					miRNAId = miRNAId + "-left+" + Math.abs(left);
				
				right = position - Integer.valueOf(miRNAIsoform[1]);
				if (right >= 0)
					miRNAId = miRNAId + "-right-" + right;
				else if (right < 0)
					miRNAId = miRNAId + "-right+" + Math.abs(right);
			}
		
			if (left == 0 && right == 0)
				miRNAId = miRNA[0] + ":" + miRNAIsoform[0];
			
			//Find out the mismatch location and create a new isoform record.
			if (mismatchNum > 0)
			{
				Integer[] locations = MiRNAUtil.getInstance().getAbsoluteMismatchLocation(MDString);
				for (int k = 0; k < locations.length; k++)
				{
					int absoluteLocation = position + locations[k];
					int relativeLocation = MiRNAUtil.getInstance().getRelativeGenomeMismatchLocation(absoluteLocation, miRNA, miRNAIsoform);
					
					if (!m_ignoreMismatchLocation && Math.abs(relativeLocation) > 8)
						return;
					
					if (relativeLocation > 0)
						miRNAId = miRNAId + ":M+" + relativeLocation;
					else
						miRNAId = miRNAId + ":M" + (relativeLocation -1);
				}
			}
			
			//Find out the "N" location and ignore the reads which location of "N" occurred within 5p.
			rawSequence = MiRNAUtil.getInstance().reverseSequence(rawSequence, miRNA);
					
			if (rawSequence.indexOf("N") > -1)
			{
				int absoluteLocation = position + rawSequence.indexOf("N");
				int relativeLocation = MiRNAUtil.getInstance().getRelativeNLocation(absoluteLocation, miRNA, miRNAIsoform);
				if (Math.abs(relativeLocation) < 8)
					return;			
			}
			
			String newSAMRecord = MiRNAUtil.getInstance().processSAMRecord(stringLine, miRNAId, miRNA, position, rawSequence);
			if (newSAMRecord == null)
				return;
			
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
				ann.add(String.valueOf(left));
				ann.add(String.valueOf(right));
				ann.add(rawSequence);
				
				m_miRNAAnnotationList.put(miRNAId, ann);
				if (!m_miRNAReadsList.containsKey(miRNAId))
				{
					m_miRNAReadsList.put(miRNAId, newSAMRecord);
					m_miRNAReadsIDList.put(miRNAId, readID);
				}
			}
			else 
			{
				ArrayList<Integer> countsList = m_miRNACountsList.get(miRNAId);
				countsList.set(fileIndex, countsList.get(fileIndex) + 1);
				
				ArrayList<String> ann = m_miRNAAnnotationList.get(miRNAId);
				int numN = MiRNAUtil.getInstance().countN(rawSequence);
				
				if (numN < MiRNAUtil.getInstance().countN(ann.get(7)))
				{
					ann.set(7, rawSequence);		
					m_miRNAReadsList.put(miRNAId, newSAMRecord);
					m_miRNAReadsIDList.put(miRNAId, readID);
				}
			}
		}
		else
		{
			m_totalIgnoredReads++;
			
			rawSequence = MiRNAUtil.getInstance().reverseSequence(rawSequence, miRNA);
			
			String newSAMRecord = MiRNAUtil.getInstance().processSAMRecord(stringLine, miRNA[0] + ":undefined", miRNA, position, rawSequence);
			if (newSAMRecord == null)
				return;
			if (!m_miRNACountsList.containsKey(miRNA[0] + ":undefined"))
			{
				ArrayList<Integer> countsList = new ArrayList<Integer>();
				for (int i = 0; i < m_fileNameList.size(); i++)
					countsList.add(0);
				
				countsList.set(fileIndex, 1);
				
				m_miRNACountsList.put(miRNA[0] + ":undefined", countsList);
				
				ArrayList<String> ann = new ArrayList<String>();
				ann.add(miRNA[0]);
				ann.add(miRNA[0]);
				ann.add(miRNA[3]);
				ann.add(miRNA[4]);
				ann.add(String.valueOf(position));
				ann.add(String.valueOf(position + length - 1));
				ann.add(String.valueOf(0));
				ann.add(String.valueOf(0));
				ann.add(rawSequence);
				m_miRNAAnnotationList.put(miRNA[0] + ":undefined", ann);
				
				if (!m_miRNAReadsList.containsKey(miRNA[0] + ":undefined"))
				{
					m_miRNAReadsList.put(miRNA[0] + ":undefined", newSAMRecord);
					m_miRNAReadsIDList.put(miRNA[0] + ":undefined", readID);
				}
			}
			else
			{
				ArrayList<Integer> countsList = m_miRNACountsList.get(miRNA[0] + ":undefined");
				countsList.set(fileIndex, countsList.get(fileIndex) + 1);	
				
				ArrayList<String> ann = m_miRNAAnnotationList.get(miRNA[0] + ":undefined");
				int numN = MiRNAUtil.getInstance().countN(rawSequence);
				
				if (numN < MiRNAUtil.getInstance().countN(ann.get(7)))
				{
					ann.set(7, rawSequence);	
					m_miRNAReadsList.put(miRNA[0] + ":undefined", newSAMRecord);
					m_miRNAReadsIDList.put(miRNA[0] + ":undefined", readID);
				}
			}			
			return;
		}
	}
	
	public void loadDesignTable(String designFile)
	{	
		try {
			m_miRNAUtil.loadDesignFile(designFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

