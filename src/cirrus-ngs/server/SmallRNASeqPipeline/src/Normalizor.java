import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;
/**
 * The class is to normalize the count table and output the normalized count table
 * @author guorong
 *
 */
public class Normalizor {
	public ArrayList<ArrayList<String>> m_miRNACounts = new ArrayList<ArrayList<String>>();
	public HashMap<String, ArrayList<Float>> m_miRNANormCounts = new HashMap<String, ArrayList<Float>>();
	public ArrayList<Float> m_totalReads = new ArrayList<Float>();
	public ArrayList<String> m_miRNAHead = new ArrayList<String>();
	public int m_annIndex = 9; // Normally the index should be 9 or 1 for TCGA
	
	public static void main(String[] args) throws Exception
	{
		String matrixFile = "/Users/Guorong/Workspace/miRNA-data/20130320_Serum_TG2_C1TNDACXX_Aaron/output_new_group/mirRNA.all.isoform_counts.txt";
		
		Normalizor normalizor = new Normalizor();
		normalizor.parseFile(matrixFile);
		normalizor.outputMatrix(matrixFile.substring(0, matrixFile.length() - 4) + ".normalized.txt");
	}
	
	public void parseFile(String file) throws Exception
	{
		try
		{
			String stringLine;
			int sampleNum = 0;
			BufferedReader in = new BufferedReader(new FileReader(file));

			while ((stringLine = in.readLine()) != null)
			{
				if (stringLine.startsWith("miRNA_name") || stringLine.startsWith("miRNA_ID"))
				{
					sampleNum = preprocessMatrixFile(stringLine);
					for (int i = 0; i < sampleNum; i++)
						m_totalReads.add(0.0f);
					continue;
				}
				
				processMatrixFile(stringLine, sampleNum);
			}

			if (in != null)
				in.close();

		} catch (Exception e)
		{
			System.err.println(e.getMessage());
		}
	}
	
	private int preprocessMatrixFile(String stringLine)
	{
		int column = 0;
		int sampleNum = 0;
		StringTokenizer st = new StringTokenizer(stringLine, "\t");
		
		while (st.hasMoreTokens())
		{
			String stringValue = st.nextToken();
			m_miRNAHead.add(stringValue);
			if (column > m_annIndex)
				sampleNum++;
			
			column++;
		}
		
		return sampleNum;
	}
	
	private void processMatrixFile(String stringLine, int sampleNum)
	{
		int column = 0;
		StringTokenizer st = new StringTokenizer(stringLine, "\t");
		ArrayList<String> isoform = new ArrayList<String>();
		
		while (st.hasMoreTokens())
		{
			String stringValue = st.nextToken();
			isoform.add(stringValue);
			
			if (column > m_annIndex)
			{
				int index = column - m_annIndex - 1;
				float count = Float.valueOf(stringValue);
				m_totalReads.set(index, m_totalReads.get(index) + count);		
			}
			column++;
		}
		
		m_miRNACounts.add(isoform);
	}
	
	public void outputMatrix(String outputFile)
	{
		try
		{
			FileWriter filewriter = new FileWriter(outputFile, true);
			for (int i = 0; i < m_miRNAHead.size(); i ++)
				filewriter.write(m_miRNAHead.get(i) + "\t");
			filewriter.write("\n");
			
			for (int i = 0; i < m_miRNACounts.size(); i++)
			{
				ArrayList<String> miRNA = m_miRNACounts.get(i);
				ArrayList<Float> normalizedCounts = new ArrayList<Float>();
				
				for (int j = 0; j < miRNA.size(); j++)
				{
					if (j < (m_annIndex + 1))
						filewriter.write(miRNA.get(j) + "\t");
					else
					{
						float count = 1000000 * Float.valueOf(miRNA.get(j))/m_totalReads.get(j - (m_annIndex + 1)) ;
						normalizedCounts.add(count);
						
						filewriter.write(count+ "\t");
					}
				}
				
				m_miRNANormCounts.put(miRNA.get(0), normalizedCounts);
				
				filewriter.write("\n");
			}
			filewriter.close();
		} catch (Exception e)
		{
			System.err.println(e.getMessage());
		}
	}
}
