import java.util.ArrayList;
/**
 * The class is to execute the S3 upload.
 * @author guorong
 *
 */
public class S3UploadExecutor extends CommandExecutor {
	public String m_separator = System.getProperty("file.separator");
	public String m_shellScript = "";
	public String m_S3UploadURL = "";
	public String m_dataFolder = "";

	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add(m_shellScript);
		commandArray.add(m_S3UploadURL);
		commandArray.add(m_dataFolder);

		return commandArray.toArray(new String[commandArray.size()]);
	}

	public void setShellFile(String fileName) {
		m_shellScript = fileName;
	}

	public void setS3UploadURL(String S3UploadURL) {
		m_S3UploadURL = S3UploadURL;
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
