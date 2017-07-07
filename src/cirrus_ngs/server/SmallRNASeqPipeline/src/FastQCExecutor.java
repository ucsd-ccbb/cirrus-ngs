import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
/**
 * The class is to execute the fastQC package.
 * @author guorong
 *
 */
public class FastQCExecutor extends CommandExecutor {
	public String m_separator = System.getProperty("file.separator");
	public String m_shellScript = "";
	public String m_disableLog = "";
	public String m_suffix = "";
	public ArrayList<String> m_sequenceFiles = new ArrayList<String>();

	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		commandArray.add("inputFile=" + filePath + fileName);

		commandArray.add(m_shellScript);

		return commandArray.toArray(new String[commandArray.size()]);
	}

	public void setShellFile(String fileName) {
		m_shellScript = fileName;
	}

	public void run(String shell, String fileName) 
	{
		setShellFile(shell);

		System.out.println("FastQC is processing:" + fileName);
		execute("", fileName);
	}
}
