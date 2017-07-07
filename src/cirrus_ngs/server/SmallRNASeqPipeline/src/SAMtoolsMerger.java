import java.util.ArrayList;
/**
 * The class is to merge small bam files into a big one bam file.
 * @author guorong
 *
 */
public class SAMtoolsMerger extends CommandExecutor {
	public String m_shellScript = "";
	public String m_fastqFile = null;
	public int m_splitNum = 10;
	
	public String[] assembleCommands(String filePath, String fileName) {

		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("qsub");
		commandArray.add("-v");
		
		String inputFiles = "inputFiles=";
		for (int i = 0; i < m_splitNum; i++)
			inputFiles = inputFiles + m_fastqFile + "." + i + ".bam ";
		
		commandArray.add("fastqFile=" + m_fastqFile + ",headSAMFile=" + m_fastqFile + ".0.sam,outputBAMFile=" + m_fastqFile + ".bam," 
				+ inputFiles.substring(0, inputFiles.length() - 1));
		commandArray.add(m_shellScript);

		return commandArray.toArray(new String[commandArray.size()]);
	}

	public void setShellFile(String fileName) {
		m_shellScript = fileName;
	}
	
	public void setFastqFile(String fastqFile) {
		m_fastqFile = fastqFile;
	}
	
	public void setSplitNum(int splitNum)
	{
		m_splitNum = splitNum;
	}
	
	public void run(String shell, String fastqFile) 
	{
		setShellFile(shell);
		setFastqFile(fastqFile);
		
		execute("", "");
	}	
}
