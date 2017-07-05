package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * MergeBAMExecutor merges splitted BAM files.
 * 
 * @author Guorong Xu
 */
public class MergeBAMExecutor extends CommandExecutor {
	
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**an arraylist of strings stores input file names*/
	public ArrayList<String> m_inputFileNames = null;
	/**a string stores the output file name*/
	public String m_outputFileName = "";
	/**a string stores the log path*/
	public String m_logFile = "";

	/**
	 * Assembles a shell command with specified shell script path, output file
	 * name, log path and input file names.
	 * The function overrides the assembleCommands from the class 
	 * CommandExecutor
	 * @param filePath		unused in the function
	 * @param fileName		unused in the function	
	 * @return String[]		a string array includes the command 
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(m_outputFileName);
		commandArray.add(m_logFile);
		for (int i = 0; i < m_inputFileNames.size(); i++)
			commandArray.add(m_inputFileNames.get(i));
		
		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * Stores the input path as the shell script path.
	 * @param fileName	The path will be stored in the field m_shellScript
	 */
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	/**
	 * Stores input file names
	 * @param inputFileNames	The Arraylist of strings will be stored in 
	 * 							m_inputFileNames
	 */
	public void setInputFileNames(ArrayList<String> inputFileNames)
	{
		m_inputFileNames = inputFileNames;
	}
	
	/**
	 * Stores the output file name
	 * @param outputFileName	The string will be stored in m_outputFileName
	 */
	public void setOutputFileName(String outputFileName)
	{
		m_outputFileName = outputFileName;
	}
	
	/**
	 * Stores the input path as the log path.
	 * @param logFile	The path will be stored in the field m_logFile
	 */
	public void setLogFile(String logFile)
	{
		m_logFile = logFile;
	}
}
