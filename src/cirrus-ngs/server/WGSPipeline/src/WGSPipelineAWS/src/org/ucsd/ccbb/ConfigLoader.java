package org.ucsd.ccbb;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

/**
 * ConfigLoader loads the configuration file which contains all paths of the 
 * tools and parameters used in the pipeline.
 * 
 * @author Guorong Xu
 */
public class ConfigLoader {
	/**a HashMap to store all the tools and corresponding paths*/
	public HashMap<String, String>  m_configsList = new HashMap<String, String>();
	
	/**
	 * Reads in a file name and calls the function processPathFile to process
	 * each line in the input file
	 * @param file  each line in the file needs to be processed
	 * @throws Exception
	 */
	public void parsePathFile(String file) throws Exception
	{
		String stringLine = null;
		try
		{		
			BufferedReader in = new BufferedReader(new FileReader(file));
			
			/*Process all lines in the file.
			 *Ignore the line started with "#" 
			 */
			while ((stringLine = in.readLine()) != null)
			{
				/*ignore comments*/
				if (stringLine.startsWith("#"))
					continue;

				processPathFile(stringLine);
			}

			if (in != null)
				in.close();
		} catch (Exception e)
		{
			System.err.println(stringLine);
		}
	}
	
	/**
	 * Stores values of tools and corresponding paths to m_configsList.
	 * It will be called by the function parsePathFile.
	 * @param stringLine	a string consists of a tool and the corresponding 
	 * 					 	path
	 * @throws IOException
	 */
	private void processPathFile(String stringLine) throws IOException
	{
		/*Construct a string tokenizer with the delimiter "\t". */
		StringTokenizer st = new StringTokenizer(stringLine, "\t");
		
		int columnth = 0;
		String tool = null;
		String path = null;
			
		/*Ignore the empty string*/
		if ((stringLine.length() == 0))
			return;

		/*Save the first token to tool.
		 *Save the second token to path.
		 */
		while (st.hasMoreTokens())
		{
			String stringValue = st.nextToken();

			switch (columnth)
			{
			case 0:
				tool = stringValue;
				break;
			case 1:
				path = stringValue;
				break;
			default:
				;
			}
			columnth++;
		}
		
		m_configsList.put(tool, path);
	}
	
	/**
	 * Aims to return mapped values associated with the key.
	 * If the input key is null or there is no such key maintained in the map,
	 * return the defaultValue
	 * @param key  the value associated with the key will be returned
	 * @param defaultValue  the value will be returned if the key is null or
	 * 					    the key is not maintained in the map
	 * @return	the value associated with the input key or the defaultValue
	 */
	public String get(String key, String defaultValue)
	{
		if (key != null && m_configsList.containsKey(key))
			return m_configsList.get(key);	
		else 
			return defaultValue;			
	}
}
