/** 
 *
 * $Author: Guorong Xu
 * $Revision: 2.7
 * $Date: 2011/09/19 
 *
 * This file is part of SAMMate.
 * Copyright (C) 2008-2011 Guorong Xu, All Rights Reserved. 
 * 
 */
import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * The class is an executor to run Linux commands.
 * @author guorong
 *
 */
public abstract class CommandExecutor {

    public String m_separator = System.getProperty("file.separator");
    private String[] m_commandLine = null;
    public String m_returnValue = "0";

    public boolean execute(String filePath, String fileName) {
        try {
            Process process = null;

            m_commandLine = assembleCommands(filePath, fileName);

            String commands = "";
            for (int i = 0; i < m_commandLine.length; i++) {
                commands = commands + m_commandLine[i] + " ";
            }

            process = Runtime.getRuntime().exec(m_commandLine);

            String line;
            BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));

            while ((line = in.readLine()) != null) {
                System.out.println(line);
            }

            in = new BufferedReader(new InputStreamReader(process.getErrorStream()));

            while ((line = in.readLine()) != null) {
                System.out.println(line);
            }

            in.close();

            int exitVal = process.waitFor();

            process.destroy();
            process = null;

            if (exitVal < 0) {
            	System.out.println("Exited with error code " + exitVal);
            	System.out.println("The command can not be successfully executed.");
                return false;
            }
            if (exitVal > 0) {
            	System.out.println("Exited with error code " + exitVal);
            	System.out.println("The command can not be successfully executed.");
                return false;
            }

        } catch (Exception e) {
        	System.out.println("The command can not be successfully executed.");
            System.out.println(e);
            return false;
        }

        return true;
    }

    public abstract String[] assembleCommands(String filePath, String fileName);
}
