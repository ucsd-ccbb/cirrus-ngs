package org.ucsd.ccbb;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

/**
 * PipelineUtil Lists all files with specified suffix in a specified directory,
 * sorts the arraylist of file directories or names in ascending order.
 * 
 * @author Guorong Xu
 */
public class PipelineUtil {

	/**an arraylist of strings stores the path of files whose name ends with 
	 * specified suffix
	 */
	public ArrayList<String> m_fileList = new ArrayList<String>();
	/**an arraylist of strings stores file names with specified suffix*/
	public ArrayList<String> m_fileNameList = new ArrayList<String>();
	
	/**
	 * Lists all files with specified suffix in a specified 
	 * directory
	 * @param folder	a directory that contains aim files
	 * @param suffix	try to find files whose names include the suffix
	 */
	public void listAllFiles(String folder, String suffix) {
		/* Returns an array of abstract pathnames denoting the files in the 
		 * directory denoted by this abstract pathname.
		 */
		File[] listOfFiles = new File(folder).listFiles();
		if (listOfFiles == null || listOfFiles.length == 0) {
			return;
		}

		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isFile()) {
				String fileName = listOfFiles[i].getName();
				if (fileName.endsWith(suffix)) {
					m_fileList.add(listOfFiles[i].getPath().substring(0, listOfFiles[i].getPath().indexOf(suffix)));
					m_fileNameList.add(listOfFiles[i].getName());
					System.out.println(listOfFiles[i].getPath().substring(0, listOfFiles[i].getPath().indexOf(suffix)));
				}
			} else if (listOfFiles[i].isDirectory()) {
				listAllFiles(listOfFiles[i].getPath(), suffix);
			}
		}
	}

	/**
	 * Sorts the arraylist of file directories in ascending order
	 * @return	an arraylist with sorted file directories
	 */
	public ArrayList<String> getFileList() {
		Collections.sort(m_fileList);
		return m_fileList;
	}
	
	/**
	 * Sorts the arraylist of file names in ascending order
	 * @return	an arraylist with sorted file names
	 */
	public ArrayList<String> getFileNameList() {
		Collections.sort(m_fileNameList);
		return m_fileNameList;
	}
}
