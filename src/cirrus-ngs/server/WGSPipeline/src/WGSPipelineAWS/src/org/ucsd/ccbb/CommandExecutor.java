package org.ucsd.ccbb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStreamReader;

/** 
 * CommandExecutor is an abstract class which works as an executor to run 
 * Linux commands.
 * BAM2FastqExecutor, S3DownLoadExecutor are extended from this class.
 * Copyright (C) 2014-2015 Guorong Xu, All Rights Reserved.
 * 
 * @author: Guorong Xu
 * @date: 2015/05/19 
 */
public abstract class CommandExecutor {
	/**get the character that separates components of a file path*/
    public String m_separator = System.getProperty("file.separator");
    /**a string array includes the command*/
    private String[] m_commandLine = null;
    /***/
    public String m_returnValue = "0";

    /**
     * The function execute the commands associated with the input file path 
     * and file name. Print out the corresponding normal output and error 
     * output. 
     * @param filePath	the file path may be used to assemble a command
     * @param fileName	the file name may be used to assemble a command
     * @return	 true if the subprocess represented by this Process object 
     * 			 terminate normally, false otherwise.  
     */
    public boolean execute(String filePath, String fileName) {
        String commands = "";
    		
        try {
            Process process = null;

            m_commandLine = assembleCommands(filePath, fileName);

            /*convert the string array to a string*/
            for (int i = 0; i < m_commandLine.length; i++) 
            {
                commands = commands + m_commandLine[i] + " ";
            }
            System.out.println("The running command: " + commands);
            /*return the runtime object associate with ?*/
            process = Runtime.getRuntime().exec(m_commandLine);

            String line;
            /*read in the input stream connected to the normal output of the subprocess*/ 
            BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));

            while ((line = in.readLine()) != null) 
            		System.out.println(line);

            /*read in the input stream connected to the error output of the subprocess*/ 
            in = new BufferedReader(new InputStreamReader(process.getErrorStream()));

            while ((line = in.readLine()) != null) 
            		System.err.println(line);

            in.close();
            /* get the exit value of the subprocess represented by this Process
             * object. 
             * By convention, the value 0 indicates normal termination*/
            int exitVal = process.waitFor();
            
            /*kill the subprocesses*/
            process.destroy();
            process = null;

            if (exitVal != 0) 
            {
	            	System.err.println("Exited with error code " + exitVal);
	            	System.err.println("The running command: " + commands);
	            	System.err.println("The command can not be successfully executed.");
                return false;
            }
      
        } catch (Exception e) {
        		System.err.println("The running command: " + commands);
        		System.err.println("The command can not be successfully executed.");
            System.err.println(e);
            return false;
        }

        return true;
    }

    /**
     * The abstract class will be overrided in subclasses.
     * @param filePath		the file path may be used to assemble a command
     * @param fileName		the file name may be used to assemble a command
     * @return String[]		a string array includes the command
     */
    public abstract String[] assembleCommands(String filePath, String fileName);
}
