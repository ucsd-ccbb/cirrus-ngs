package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * DeleteWorkspaceExecutor deletes workspace.
 * The class extends the class CommandExecutor.
 * 
 * @author Guorong Xu
 */
public class DeleteWorkspaceExecutor extends CommandExecutor {
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the input file path*/
	public String m_inputFile = "";
	
	/**
	 * Assembles a shell command with shell script path and input file path
	 * The function overrides the function assembleCommands() from class 
	 * CommandExecutor.
	 * @param filePath	unused in the function
	 * @param fileName	unused in the function
	 * @return a string array includes the command 
	 */
	public String[] assembleCommands(String filePath, String fileName) {
		ArrayList<String> commandArray = new ArrayList<String>();

		commandArray.add("sh");
		commandArray.add(m_shellScript);
		commandArray.add(m_inputFile);
		
		return commandArray.toArray(new String[commandArray.size()]);
	}
	
	/**
	 * This function set the input path as the shell script path.
	 * @param fileName	The path will be set to the field m-shellScript
	 */
	public void setShellFile(String fileName)
	{
		m_shellScript = fileName;
	}
	
	/**
	 * This function set the file path as the input file path
	 * @param inputFile		The file path will be set to be the input file path
	 */
	public void setInputFile(String inputFile)
	{
		m_inputFile = inputFile;
	}
}
