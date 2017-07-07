import java.util.ArrayList;
/**
 * The class is to clean all temporary files after processing is over.
 * @author guorong
 *
 */
public class FileCleaner extends CommandExecutor {
	public String m_shellScript = "";
	public String m_fileName = "";
	public String m_outputSAMFile = "";
	
	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		commandArray.add("fileName=" + m_fileName + ",outputSAMFile=" + m_outputSAMFile);
		commandArray.add(m_shellScript);

		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	public void setFileName(String fileName)
	{
		m_fileName = fileName;
	}
	
	public void setOutputSAMFile(String outputSAMFile)
	{
		m_outputSAMFile = outputSAMFile;
	}
	
	public void run(String shellFile, String fileName, String outputSAMFile) 
	{		
		setShellFile(shellFile);
		setFileName(fileName);
		setOutputSAMFile(outputSAMFile);
		execute("", "");
	}
}
