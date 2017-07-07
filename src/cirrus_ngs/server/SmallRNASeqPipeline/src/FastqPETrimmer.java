import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
/**
 * The class is to trim the pair-end fastq files.
 * @author guorong
 *
 */
public class FastqPETrimmer extends CommandExecutor {
	
	public String m_shellScript = "";
	public ArrayList<String> m_fastq1 = new ArrayList<String>();
	public ArrayList<String> m_fastq2 = new ArrayList<String>();
	public String m_fastqFile1 = "";
	public String m_fastqFile2 = "";
	public String m_suffix = "";
	
	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		commandArray.add("type=PE," + "inputFile1=" + m_fastqFile1+",inputFile2=" + m_fastqFile2);
		
		commandArray.add(m_shellScript);

		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	public void setSuffix(String suffix)
	{
		m_suffix = suffix;
	}
}
