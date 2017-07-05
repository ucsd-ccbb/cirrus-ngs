import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
/**
 * The class is to merge the novoalign log files.
 * @author guorong
 *
 */
public class MergeLogsExecutor extends CommandExecutor {
	public String m_separator = System.getProperty("file.separator");
	public String m_shellScript = "";

	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		commandArray.add("inputFolder=" + filePath + fileName + ",outputFile=" + filePath + fileName 
				+ System.getProperty("file.separator") + "output" + System.getProperty("file.separator") + "alignment.log");

		commandArray.add(m_shellScript);

		return commandArray.toArray(new String[commandArray.size()]);
	}

	public void setShellFile(String shellScript) {
		m_shellScript = shellScript;
	}

	public void run(String shell, String filePath) 
	{
		setShellFile(shell);

		System.out.println("MergeLogs is processing:" + filePath);
		execute("", filePath);
	}
}
