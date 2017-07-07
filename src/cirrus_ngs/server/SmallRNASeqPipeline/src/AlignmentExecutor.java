import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
/**
 * The class is the executor for alignment by novoalign.
 * @author guorong
 *
 */
public class AlignmentExecutor extends CommandExecutor {
	
	public String m_shellScript = "";
	public ArrayList<String> m_fastq = new ArrayList<String>();
	public String m_fastqFile = "";
	public String m_logFilePath = "";
	public String m_suffix = "";
	public ConfigLoader m_configLoader = null;
	
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		
		commandArray.add("-o");
		commandArray.add(m_logFilePath);
		commandArray.add("-e");
		commandArray.add(m_logFilePath);
		
		commandArray.add("-v");
		commandArray.add("fastqFile=" + m_fastqFile + ",outputFile=" + fileName);
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
	
	public void setConfigLoader(ConfigLoader configLoader)
	{
		m_configLoader = configLoader;
	}
	
	public void setLogFilePath(String logFilePath)
	{
		m_logFilePath = logFilePath;
	}
	
	public void runBowtie(String shellFile, String inputFile, String fastqFileName) throws InterruptedException {	
		
		setShellFile(shellFile);
		setLogFilePath(inputFile);
	
		String outputFileName = inputFile;		
		m_fastqFile = fastqFileName;
		
		System.out.println("Bowtie is processing:" + fastqFileName);
		execute("", outputFileName);
	}
}
