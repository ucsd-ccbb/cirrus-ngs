import java.util.ArrayList;
/**
 * The class is to execute the S3 download.
 * @author guorong
 *
 */
public class S3DownloadExecutor extends CommandExecutor {
	public String m_separator = System.getProperty("file.separator");
	public String m_shellScript = "";
	public String m_S3DownloadURL = "";
	public String m_dataFolder = "";

	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add(m_shellScript);
		commandArray.add(filePath + fileName);
		commandArray.add(m_S3DownloadURL);
		commandArray.add(m_dataFolder);

		return commandArray.toArray(new String[commandArray.size()]);
	}

	public void setShellFile(String fileName) {
		m_shellScript = fileName;
	}

	public void setS3DownloadURL(String S3DownloadURL) {
		m_S3DownloadURL = S3DownloadURL;
	}
	
	public void setDataFolder(String dataFolder) {
		m_dataFolder = dataFolder;
	}
	
	public void run(String shell, String fileName) 
	{
		setShellFile(shell);
		execute("", fileName);
	}
}
