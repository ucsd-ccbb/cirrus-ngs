import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;
/**
 * The class is sort a matrix table according to multiple conditions.
 * @author guorong
 *
 */
class MatrixSorter {
	public void printAL(Object[] fields, String outputFileName) {
		File outputFile = new File(outputFileName + ".sort");
		FileWriter filewriter;
		try {
			filewriter = new FileWriter(outputFile, true);
			for (int i = 0; i < fields.length; i++) {
				StringTokenizer stk = new StringTokenizer(String.valueOf(fields[i]) , ",");
				while (stk.hasMoreTokens())
				{
					String stringValue = stk.nextToken();
					filewriter.write(stringValue + "\t");
				}
				filewriter.write("\n");		
			}
			filewriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void sort(String inputFileName) {
		ArrayList<Field> al = new ArrayList<Field>();
		try {
			BufferedReader br = new BufferedReader(
					new FileReader(inputFileName));
			String line = br.readLine();
			while (line != null) {
				if (line.startsWith("#"))
				{
					line = br.readLine();
					continue;
				}

				StringTokenizer stk = new StringTokenizer(line, "\t");
				al.add(new Field(stk.nextToken().trim(),
						stk.nextToken().trim(), stk.nextToken().trim(), stk.nextToken().trim(), stk.nextToken().trim(),
						stk.nextToken().trim()));
				line = br.readLine();
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Object[] fields = al.toArray();
		Arrays.sort(fields, new OrderByComparator());

		printAL(fields, inputFileName);
	}
	
	public static void main(String[] args) throws Exception {	
		MatrixSorter sorter = new MatrixSorter();
		
		sorter.sort("/Users/Guorong/Workspace/miRNA-data/20130320_PBMC_TG2_C1TNDACXX_Aaron/mirRNA.all.isoform_counts_total.txt");

	}
}