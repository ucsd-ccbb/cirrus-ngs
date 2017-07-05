import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
/**
 * The class is to call rnapipeline.sh script which is to convert the artificial chromosome name to regular chromosome name,
 * to sort the SAM file and to convert SAM file to BAM file.
 * @author guorong
 *
 */
public class FormatConversion extends CommandExecutor {
	
	public String m_shellScript = "";
	public ArrayList<String> m_samFiles = new ArrayList<String>();
	public String m_samSuffix = "";
	public String m_outputFolder = "";
	
	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		commandArray.add("inputFile=" + filePath + fileName);
		commandArray.add(m_shellScript);

		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	public void setSAMSuffix(String suffix)
	{
		m_samSuffix = suffix;
	}

	public void run(String shell, String inputFile) 
	{
		setShellFile(shell);
		
		System.out.println("FormatConversion is processing:" + inputFile + ".sam");
		execute("", inputFile);
	}	
}
