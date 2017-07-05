package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * GroupVCFExecutor groups individual VCF files into one family.
 * 
 * @author Guorong Xu
 */
public class GroupVCFExecutor extends CommandExecutor {
	
	/**a string stores the shell script path*/
	public String m_shellScript = "";
	/**a string stores the group name*/
	public String m_groupName = "";
	/**an arraylist stores a list of VCF files*/
	public ArrayList<String> m_vcfList = new ArrayList<String>();
	/**a string stores the log path*/
	public String m_logFile = "";
	
	/**
	 * Assembles a shell command with specified shell script path, group name, 
	 * log path and a list of VCF files.
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
		commandArray.add(m_groupName);
		commandArray.add(m_logFile);
		for (int i = 0; i < m_vcfList.size(); i ++)
			commandArray.add(m_vcfList.get(i));

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
	 * Stores the input string as the group name
	 * @param groupName		The string will be stored in the field m_groupName
	 */
	public void setGroupName(String groupName)
	{
		m_groupName = groupName;
	}
	
	/**
	 * Stores a list of VCF files
	 * @param vcfList	a list of VCF files will be stored in the field 
	 * 					m_vcfList
	 */
	public void setVcfList(ArrayList<String> vcfList)
	{
		m_vcfList = vcfList;
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
