import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
/**
 * The class is to trim the single-end fastq files.
 * @author guorong
 *
 */
public class FastqSETrimmer extends CommandExecutor {
	
	public String m_shellScript = "";
	public ArrayList<String> m_fastq = new ArrayList<String>();
	public String m_fastqFile = "";
	public String m_suffix = "";
	
	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		commandArray.add("type=SE," + "inputFile=" + m_fastqFile);
		
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
