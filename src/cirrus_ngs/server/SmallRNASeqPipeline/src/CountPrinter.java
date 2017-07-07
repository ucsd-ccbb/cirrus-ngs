import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

/**
 * The class is to print out all count table
 * @author guorong
 *
 */
public class CountPrinter {
	public ArrayList<String[]> m_totalCountsRecord = null;
	public MiRNAUtil m_miRNAUtil = MiRNAUtil.getInstance();
	public Normalizor m_normalizor = null;
	
	public void execute(String inputPath) throws Exception
	{
		String label = "all";
		String countsType = "average";

		String countFile = inputPath + "mirRNA." + label + ".isoform_counts_" + countsType + ".txt.sort";
		String matrixFile = inputPath + "mirRNA." + label + ".isoform_counts.txt";
		String groupFile = inputPath + "group.txt";
		
		Normalizor normalizor = new Normalizor();
		normalizor.parseFile(matrixFile);
		normalizor.outputMatrix(matrixFile.substring(0, matrixFile.length() - 4) + ".normalized.txt");
		
		CountPrinter countPrinter = new CountPrinter();
		countPrinter.loadDesignTable(groupFile);
		countPrinter.setNormalizor(normalizor);
		countPrinter.loadCountsFile(countFile);
		countPrinter.outputNormalizedCounts(countFile.substring(0, countFile.length() - 9) + ".normalized.txt", true);
		
		countsType = "total";
		countFile = inputPath + "mirRNA." + label + ".isoform_counts_" + countsType + ".txt.sort";
		
		countPrinter = new CountPrinter();
		countPrinter.loadDesignTable(groupFile);
		countPrinter.setNormalizor(normalizor);
		countPrinter.loadCountsFile(countFile);
		countPrinter.outputNormalizedCounts(countFile.substring(0, countFile.length() - 9) + ".normalized.txt", false);
	}
	
	/**
	 * The method is to load the design table (group table).
	 * @param designFile
	 */
	public void loadDesignTable(String designFile)
	{	
		try {
			m_miRNAUtil.loadDesignFile(designFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void setNormalizor(Normalizor normalizor)
	{
		m_normalizor = normalizor;
	}
	
	/**
	 * The method is to load the counts files including the annotation fields
	 * @param countsFile
	 */
	public void loadCountsFile(String countsFile)
	{
		try
		{
			m_totalCountsRecord = new ArrayList<String[]>();
			
			String stringLine;
			BufferedReader in = new BufferedReader(new FileReader(countsFile));

			while ((stringLine = in.readLine()) != null)
			{
				StringTokenizer st = new StringTokenizer(stringLine, "\t");

				int columnth = 0;
				String fileName = null;
				String familyName = null;
				String isoName = null;
				String readID = null;
				String type = null;
				String number = null;
				
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
						number = stringValue;
						break;
					default:
						;
					}
					columnth++;
				}
				
				m_totalCountsRecord.add(new String[]{fileName, familyName, isoName, readID, type, number});
			}
			
			if (in != null)
				in.close();
		} catch (Exception e)
		{
		}
	}
	
	/**
	 * The method is to output normalized counts table. 
	 * The count value should be followed after the annotation matrix.
	 * @param outputPath
	 * @param isAverage
	 */
	public void outputNormalizedCounts(String outputPath, boolean isAverage)
	{
		try
		{
			File outputFile = new File(outputPath);
			FileWriter filewriter = new FileWriter(outputFile, true);
			
			HashMap<String, ArrayList<Float>> miRNANormCounts = m_normalizor.m_miRNANormCounts;
			ArrayList<String> miRNAHead = m_normalizor.m_miRNAHead;
			HashMap<String, String> groupList = m_miRNAUtil.m_groupList;
			HashMap<String, ArrayList<Integer>> groupIndex = new HashMap<String, ArrayList<Integer>>();
			
			for (int i = (m_normalizor.m_annIndex + 1); i < miRNAHead.size(); i++)
			{
				if (groupList.containsKey(miRNAHead.get(i)))
				{
					String type = groupList.get(miRNAHead.get(i));
					if (groupIndex.containsKey(type))
					{
						ArrayList<Integer> indexList = groupIndex.get(type);
						indexList.add(i - (m_normalizor.m_annIndex + 1));						
					}
					else
					{
						ArrayList<Integer> indexList = new ArrayList<Integer>();
						indexList.add(i - (m_normalizor.m_annIndex + 1));	
						groupIndex.put(type, indexList);
					}
				}
			}
			
			for (int i = 0; i < m_totalCountsRecord.size(); i++)
			{
				String[] record = m_totalCountsRecord.get(i);
				
				// fileName/isoID, familyName, iosName, readID, type, number
				if (miRNANormCounts.containsKey(record[0]))
				{					
					String type = record[4];
					ArrayList<Float> counts = miRNANormCounts.get(record[0]);
					if (groupIndex.containsKey(type))
					{
						ArrayList<Integer> indexList = groupIndex.get(type);
						float normalizedCount = 0;
						for (int j = 0; j < indexList.size(); j++)
							normalizedCount = normalizedCount + counts.get(indexList.get(j));
						
						if (isAverage)
							record[5] = String.valueOf(normalizedCount/indexList.size());
						else
							record[5] = String.valueOf(normalizedCount);						
					}
				}
				
				filewriter.write(record[0] + "\t" + record[1] + "\t" + record[2] + "\t" + record[3] + "\t" + record[4] + "\t" + record[5] + "\n");	
			}
			
		filewriter.close();
		
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
}
