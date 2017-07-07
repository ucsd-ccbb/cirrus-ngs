import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.regex.*;

public class MiRNAUtil {
	private static MiRNAUtil m_instance = null;
	
	public HashMap<String, String> m_groupList = new HashMap<String, String>();
	public HashMap<String, String> m_groupTypes = new HashMap<String, String>();

	private MiRNAUtil()
	{

	}

	public synchronized static MiRNAUtil getInstance()
	{
		if (m_instance == null)
		{
			m_instance = new MiRNAUtil();
		}
		return m_instance;
	}
	
	public HashMap<String, String> loadDesignFile(String file) throws Exception
	{
		try
		{
			String stringLine;
			BufferedReader in = new BufferedReader(new FileReader(file));

			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("#"))
					continue;
				processFile(stringLine);
			}
			
			if (in != null)
				in.close();
		} catch (Exception e)
		{
		}
		
		return m_groupList;
	}
	
	private void processFile(String stringLine)
	{
		StringTokenizer st = new StringTokenizer(stringLine, "\t");

		int columnth = 0;
		String fileName = null;
		String type = null;
		
		if ((stringLine.length() == 0))
			return;

		while (st.hasMoreTokens())
		{
			String stringValue = st.nextToken();

			switch (columnth)
			{
			case 0:
				fileName = stringValue;
				break;
			case 1:
				type = stringValue;
				break;			
			default:
				;
			}
			columnth++;
		}
		
		m_groupList.put(fileName, type);
		
		if (type != null && !m_groupTypes.containsKey(type))
			m_groupTypes.put(type, type);
	}

	public void setGroupTypes(HashMap<String, String> groupTypes)
	{
		m_groupTypes = groupTypes;
	}

	public HashMap<String, String> getGroupTypes()
	{
		return m_groupTypes;
	}
	
	public void setGroupList(HashMap<String, String> groupList )
	{		
		m_groupList = groupList;
	}
	
	public HashMap<String, String> getGroupList()
	{
		return m_groupList;
	}
	
	public int getMismatchNum(String MDString)
	{
		if (MDString == null)
			return 0;
		
		if (MDString.startsWith("MD:"))
		{
			String[] acgt = new String[]{"A", "C", "G", "T"};
			int mismatchNum = 0;
			
			for (int i = 0; i < acgt.length; i++)
			{
				int index = MDString.indexOf(acgt[i]);
				if (index > -1)
				{
					int secondIndex = MDString.substring(index + 1).indexOf(acgt[i]);
					if (secondIndex > -1)
						mismatchNum = mismatchNum + 2;
					else
						mismatchNum = mismatchNum + 1;
				}
			}
			
			return mismatchNum;
		}
		else if (MDString.startsWith("NM:"))
		{
			int mismatchNum = Integer.valueOf(MDString.substring(MDString.length() - 1));
			return mismatchNum;
		}
		
		return 0;
	}
	
	
	public Integer[] getAbsoluteMismatchLocation(String MDString)
	{
		String[] acgt = new String[]{"A", "C", "G", "T"};
		ArrayList<Integer> locations = new ArrayList<Integer>();
		String mismatchString = MDString.substring(5, MDString.length());
		
		for (int i = 0; i < acgt.length; i++)
		{
			int index = mismatchString.indexOf(acgt[i]);
			if (index > -1)
			{
				String firstPart = mismatchString.substring(0, index);
				String secondPart = mismatchString.substring(index + 1, mismatchString.length());
				boolean hasAdded = false;
				
				for (int j = 0; j < acgt.length; j++)
				{
					int firstIndex = firstPart.indexOf(acgt[j]);
					if (firstIndex > -1)
					{
						if ((firstIndex - 2) > -1 && Character.isDigit(firstPart.charAt(firstIndex - 2)))						
							locations.add(Integer.valueOf(firstPart.substring(firstIndex - 2, firstIndex)));
						else
							locations.add(Integer.valueOf(firstPart.substring(firstIndex - 1, firstIndex)));
						
						hasAdded = true;
					}
					
					int secondIndex = secondPart.indexOf(acgt[j]);
					if (secondIndex > -1)
					{
						if ((secondIndex - 2) > -1 && Character.isDigit(secondPart.charAt(secondIndex - 2)))						
							locations.add(Integer.valueOf(secondPart.substring(secondIndex - 2, secondIndex)));
						else
							locations.add(Integer.valueOf(secondPart.substring(secondIndex - 1, secondIndex)));
						
						hasAdded = true;
					}
				}
				
				if (!hasAdded)
				{
					if (firstPart.length() == 1 && Character.isDigit(firstPart.charAt(0)))
						locations.add(Integer.valueOf(firstPart));
					else if (firstPart.length() == 2 && Character.isDigit(firstPart.charAt(0)))
						locations.add(Integer.valueOf(firstPart));
				}
			}
		}
		
		return locations.toArray(new Integer[locations.size()]);
	}
	
	public int getRelativeMirBaseMismatchLocation(int absoluteLocation, String[] miRNA, String[] miRNAIsoform)
	{
		int relativeLocation = 0;
		
		if ("+".equalsIgnoreCase(miRNA[4]))
		{
			relativeLocation = Integer.valueOf(miRNA[1]) + absoluteLocation - Integer.valueOf(miRNAIsoform[1]) ;
			return relativeLocation;
		}
		else
		{
			relativeLocation = Integer.valueOf(miRNA[2]) - absoluteLocation - Integer.valueOf(miRNAIsoform[2]);
			return -relativeLocation;
		}
	}
	
	public int getRelativeGenomeMismatchLocation(int absoluteLocation, String[] miRNA, String[] miRNAIsoform)
	{
		int relativeLocation = 0;
		
		if ("+".equalsIgnoreCase(miRNA[4]))
		{
			relativeLocation = absoluteLocation - Integer.valueOf(miRNAIsoform[1]) + 1;
			return relativeLocation;
		}
		else
		{
			relativeLocation = absoluteLocation - Integer.valueOf(miRNAIsoform[2]) - 1;
			return -relativeLocation;
		}
	}
	
	public int getRelativeNLocation(int absoluteLocation, String[] miRNA, String[] miRNAIsoform)
	{
		int relativeLocation = 0;
		
		relativeLocation = absoluteLocation - Integer.valueOf(miRNAIsoform[1]) - 1;
		return relativeLocation;
	}
	
	public int countN(String rawSequence)
	{
		Matcher m = Pattern.compile("N").matcher(rawSequence);
	    int count;
	    for (count = 0; m.find(); count++);
	    
	    return count;
	}
	
	public String reverseSequence(String rawSequence, String[] miRNA)
	{
		if ("-".equalsIgnoreCase(miRNA[4]))
		{
			String reversedSequence = "";
			for (int i = 0; i < rawSequence.length(); i++)
			{
				if (rawSequence.substring(i, i+1).equalsIgnoreCase("A"))
					reversedSequence =  "T" + reversedSequence;
				else if (rawSequence.substring(i, i+1).equalsIgnoreCase("T"))
					reversedSequence =  "A" + reversedSequence;
				else if (rawSequence.substring(i, i+1).equalsIgnoreCase("C"))
					reversedSequence =  "G" + reversedSequence;
				else if (rawSequence.substring(i, i+1).equalsIgnoreCase("G"))
					reversedSequence =  "C" + reversedSequence;
				else if (rawSequence.substring(i, i+1).equalsIgnoreCase("N"))
					reversedSequence =  "N" + reversedSequence;			
			}
			
			return reversedSequence;
		}
		
		return rawSequence;
	}
	
	public String processSAMRecord(String stringLine, String newReadID, String[] miRNA, int newPosition, String newSequence)
	{
		StringTokenizer st = new StringTokenizer(stringLine, "\t");

		int columnth = 0;
		String chrom = null;
		String flag = null;
		int position = 0;
		String rawSequence = null;
		String readID = null;
			
		if ((stringLine.length() == 0))
			return stringLine;

		while (st.hasMoreTokens())
		{
			String stringValue = st.nextToken();

			switch (columnth)
			{
			case 0:
				readID = stringValue;
				break;
			case 1:
				flag = stringValue;
				if (!"0".equalsIgnoreCase(flag) && !"16".equalsIgnoreCase(flag))
					return stringLine;
				break;			
			case 2:
				chrom = stringValue;
				if ("*".equalsIgnoreCase(chrom))
					return stringLine;
				break;
			case 3:
				position = Integer.valueOf(stringValue);
				break;
			case 9:
				rawSequence = stringValue;
				break;
			default:
				;
			}
			columnth++;
		}
		String newSAMRecord = stringLine;
		//newSAMRecord = newSAMRecord.replaceAll(readID, newReadID);
		newSAMRecord = newSAMRecord.replaceAll(chrom, miRNA[0]);
		if ("+".equalsIgnoreCase(miRNA[4]))
		{
			int startPoint = newPosition - Integer.valueOf(miRNA[1]) + 1;
			if (startPoint <= 0)
				return null;
			newSAMRecord = newSAMRecord.replaceAll(String.valueOf(position), String.valueOf(startPoint));
		}
		else
		{
			int startPoint = Integer.valueOf(miRNA[2]) - newPosition - rawSequence.length() + 2;
			if (startPoint <= 0)
				return null;
			newSAMRecord = newSAMRecord.replaceAll(String.valueOf(position), String.valueOf(startPoint));
		}
		newSAMRecord = newSAMRecord.replaceAll(rawSequence, newSequence);
		
		return newSAMRecord;
	}
}
