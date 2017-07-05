import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

/**
 * The class is to output the count table at gene level and isoform level.
 * @author guorong
 *
 */
public class MiRNAPrinter {
	public ArrayList<String> m_fileNameList = null;
	public HashMap<String, ArrayList<Integer>> m_miRNACountsList = null;
	public HashMap<String, ArrayList<String>> m_miRNAAnnotationList = null;
	public HashMap<String, ArrayList<Integer>> m_groupIndexList = new HashMap<String, ArrayList<Integer>>();
	public HashMap<String, String> m_mirRNAReadsList = new HashMap<String, String>();
	public HashMap<String, String> m_miRNAReadsIDList = new HashMap<String, String>();
	
	public MiRNAPrinter(ArrayList<String> fileNameList, HashMap<String, ArrayList<Integer>> mirRNACountsList,
			HashMap<String, ArrayList<String>> mirRNAAnnotationList, HashMap<String, String> mirRNAReadsList,
			HashMap<String, String> mirRNAReadsIDList)
	{
		m_fileNameList = fileNameList;
		m_miRNACountsList = mirRNACountsList;
		m_miRNAAnnotationList = mirRNAAnnotationList;
		m_mirRNAReadsList = mirRNAReadsList;
		m_miRNAReadsIDList = mirRNAReadsIDList;
	}
	
	public void outputIsoformCountFile(String outputPath)
	{
		File outputFile = new File(outputPath);
		FileWriter filewriter;
		try
		{
			filewriter = new FileWriter(outputFile, true);

			// Output the head of SAM file;
			filewriter.write("miRNA_name\tmiRNA_family\tmiRNA_isoform\tchrom\tstrand\tstart_position\tend_position\t" +
					"left_shifted_base\tright_shifted_base\traw_sequence" );
			for (int fileNum = 0; fileNum < m_fileNameList.size(); fileNum++)
				filewriter.write("\t" + m_fileNameList.get(fileNum));
			filewriter.write("\n");

			// Output the body of expression matrix file
			Set<String> keys = m_miRNACountsList.keySet();
			
			// Iterate over the keys
			for (String key : keys)
			{
				filewriter.write(key + "\t");
				
				ArrayList<String> ann = m_miRNAAnnotationList.get(key);
				for (int i = 0; i < ann.size(); i++)
					filewriter.write(ann.get(i) + "\t");
				
				ArrayList<Integer> countExpr = m_miRNACountsList.get(key);
				for (int i = 0; i < countExpr.size(); i++)
					filewriter.write(countExpr.get(i) + "\t");

				filewriter.write("\n");
			}

			filewriter.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public void outputMiRNACountFile(String outputPath)
	{
		createIndex();
		
		MiRNAUtil util = MiRNAUtil.getInstance();
		
		File outputFile = new File(outputPath);
		FileWriter filewriter;
		try
		{
			filewriter = new FileWriter(outputFile, true);

			// Output the head of count file;
			filewriter.write("miRNA_name\tmiRNA_family\tmiRNA_isoform\tchrom\tstrand\tstart_position\tend_position\tleft_shifted_base\t" +
					"right_shifted_base\traw_sequence" );
			
			Set<String> keys = util.m_groupTypes.keySet();
			for (String key : keys) {
				String type = util.m_groupTypes.get(key) ;
				filewriter.write("\t" + type + "_total_counts\t" + type + "_average_counts");
			}
				
			filewriter.write("\n");

			// Output the body of expression matrix file
			keys = m_miRNACountsList.keySet();
			
			// Iterate over the keys
			for (String key : keys)
			{
				filewriter.write(key + "\t");
				
				ArrayList<String> ann = m_miRNAAnnotationList.get(key);
				for (int i = 0; i < ann.size(); i++)
					filewriter.write(ann.get(i) + "\t");
				
				ArrayList<Integer> countExpr = m_miRNACountsList.get(key);
				
				Set<String> typeKeys = util.m_groupTypes.keySet();
				
				for (String typeKey : typeKeys) {
					float total = 0;
					String type = util.m_groupTypes.get(typeKey);	
					ArrayList<Integer> indexList = m_groupIndexList.get(type);
					
					for (int index = 0; index < indexList.size(); index++)
						total = total + countExpr.get(indexList.get(index));
					
					float average = total/indexList.size();
					DecimalFormat fnum = new DecimalFormat("##0.0");  
					filewriter.write(total + "\t" + fnum.format(average) + "\t");
				}

				filewriter.write("\n");
			}

			filewriter.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public void outputMiRNACountTable(String outputPath, boolean isTotal)
	{
		MiRNAUtil util = MiRNAUtil.getInstance();
		
		File outputFile = new File(outputPath);
		FileWriter filewriter;
		try
		{
			filewriter = new FileWriter(outputFile, true);

			// Output the head of count file;
			if (isTotal)
				filewriter.write("#miRNA_name\tmiRNA_family\tmiRNA_isoform\tReadID\tType\tTotal_num\n");
			else
				filewriter.write("#miRNA_name\tmiRNA_family\tmiRNA_isoform\tReadID\tType\tAverage_num\n");

			// Output the body of expression matrix file
			Set<String> keys = m_miRNACountsList.keySet();
			// Iterate over the keys
			for (String key : keys)
			{
				Set<String> typeKeys = util.m_groupTypes.keySet();
				ArrayList<String> ann = m_miRNAAnnotationList.get(key);
				ArrayList<Integer> countExpr = m_miRNACountsList.get(key);
				
				for (String typeKey : typeKeys) 
				{
					float total = 0;			
					String type = util.m_groupTypes.get(typeKey);					
					ArrayList<Integer> indexList = m_groupIndexList.get(type);
					
					filewriter.write(key + "\t" + ann.get(0) + "\t" + ann.get(1) + "\t" + m_miRNAReadsIDList.get(key) + "\t" + typeKey + "\t");	
					
					for (int index = 0; index < indexList.size(); index++)
						total = total + countExpr.get(indexList.get(index));
					
					if (!isTotal)
					{
						float average = total/indexList.size();
						DecimalFormat fnum = new DecimalFormat("##0.0");  
						filewriter.write( fnum.format(average) + "\n");
					}		
					else
						filewriter.write( total + "\n");
				}
			}

			filewriter.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	public void outputReadFile(String outputPath)
	{
		File outputFile = new File(outputPath);
		FileWriter filewriter;
		try
		{
			filewriter = new FileWriter(outputFile, true);

			// Output the body of expression matrix file
			Set<String> keys = m_mirRNAReadsList.keySet();
			// Iterate over the keys
			for (String key : keys)
				filewriter.write(m_mirRNAReadsList.get(key) + "\n");

			filewriter.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	private void createIndex()
	{
		MiRNAUtil util = MiRNAUtil.getInstance();
		HashMap<String, String> groupList = util.getGroupList();
		
		for (int i = 0; i < m_fileNameList.size(); i++)
		{
			String fileName = m_fileNameList.get(i);
			String type = groupList.get(fileName.substring(0, fileName.length() - 4));
			
			if (!m_groupIndexList.containsKey(type))
			{
				ArrayList<Integer> indexList = new ArrayList<Integer>();
				indexList.add(i);
				m_groupIndexList.put(type, indexList);
			}
			else
			{
				ArrayList<Integer> indexList = m_groupIndexList.get(type);
				indexList.add(i);
			}			
		}
	}
}
