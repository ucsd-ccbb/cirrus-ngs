import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
/**
 * The class is to split a big .gz file into many small .gz files, and the file size and file number can be set.
 * @author guorong
 *
 */
public class FileSpliter {
	private String m_fileName = "";
	private long m_fileSize = 0l;

	public void setFileName(String fileName) {
		m_fileName = fileName;
	}

	public void setFileSize(long fileSize) {
		m_fileSize = fileSize;
	}

	public long getSplittedFileSize(int fileNum) {
		long fileSize = 0;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader( new GZIPInputStream(new FileInputStream(m_fileName))));
			long lines = 0;
			while ((in.readLine()) != null)
				lines++;

			fileSize = lines / fileNum;

		} catch (Exception e) {
			System.out.println(e);
		}

		return (fileSize / 100 + 1) * 100;
	}

	public int split() {
		int fileNum = 0;
		try {
			String stringLine;
			
			BufferedReader in = new BufferedReader(new InputStreamReader( new GZIPInputStream(new FileInputStream(m_fileName))));

			long lines = 0;
			while ((stringLine = in.readLine()) != null)
				lines++;

			BufferedWriter bufferedWriter = null;
			fileNum = (int) ((lines / m_fileSize) + 1);
			File[] files = new File[fileNum];

			for (int i = 0; i < fileNum; i++) {
				files[i] = new File(m_fileName + "." + i);
				
				bufferedWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(files[i]))));
				
				in = new BufferedReader(new InputStreamReader( new GZIPInputStream(new FileInputStream(m_fileName))));
				long count = 0;
				while ((stringLine = in.readLine()) != null) {
					count++;
					if (count > m_fileSize * i && count <= m_fileSize * (i + 1))
					{	
		                bufferedWriter.write(stringLine);
		                bufferedWriter.newLine();
					}
				}

				bufferedWriter.close();
			}

			if (in != null)
				in.close();
		} catch (Exception e) {
			System.out.println(e);
		}

		return fileNum;
	}

	public static void main(String[] args) throws IOException {
		FileSpliter st = new FileSpliter();
		st.setFileName("/Users/Guorong/Workspace/miRNA-data/20130816_Tech_Dev_H0Y8DADXX_Aaron/Sample_JC_Qiagen_SpikeIn_R1.trim.fastq.gz");

		long fileSize = st.getSplittedFileSize(4);
		st.setFileSize(fileSize);
		st.split();
		System.out.println("fileSize: " + (fileSize / 100 + 1) * 100);
	}
}