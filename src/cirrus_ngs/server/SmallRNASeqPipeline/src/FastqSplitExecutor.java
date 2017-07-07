import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
/**
 * The class is to split the fastq files into small trunks.
 * @author guorong
 *
 */
public class FastqSplitExecutor extends CommandExecutor {
	public String m_separator = System.getProperty("file.separator");
	public String m_shellScript = "";
	public String m_disableLog = "";
	public String m_suffix = "";
	public int m_splitNum = 10;
	public ArrayList<String> m_sequenceFiles = new ArrayList<String>();

	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		commandArray.add("inputFile=" + filePath + fileName + ",splitNum=" + m_splitNum);

		commandArray.add(m_shellScript);

		return commandArray.toArray(new String[commandArray.size()]);
	}

	public void setShellFile(String fileName) {
		m_shellScript = fileName;
	}
	
	public void setSplitNum(int splitNum)
	{
		m_splitNum = splitNum;
	}
	
	public void setDisableLog(String disableLog) {
		m_disableLog = disableLog;
	}

	public void run(String shell, String fileName) 
	{
		setShellFile(shell);
		setDisableLog("true");

		System.out.println("FastqSplit is processing:" + fileName);
		execute("", fileName);
	}
}
