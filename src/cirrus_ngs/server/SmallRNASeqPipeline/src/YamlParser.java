

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import org.yaml.snakeyaml.Yaml;

/**
 * YamlParser is the class that parses the customer's yaml file and provides
 * informations of demands of customers.
 * A YamlParse object contains the information needed for the various callers 
 * of analysis. The information includes:
 * <ul>
 * <li>The name of the project
 * <li>Methods of analysis
 * <li>The date of the project
 * <li>The S3 upload URL
 * <li>Sample groups
 * <li>all fastq files in examples
 * <li>all BAM files in examples
 * </ul>
 * @author	Guorong Xu
 */
public class YamlParser {
	/**a string stores the name of the project*/
	public String m_projectName = "";
	/**a string stores methods of analysis"*/
	public String m_analysis = "";
	/**a string stores the project date*/
	public String m_projectDate = "";
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a HashMap stores information of all sample groups*/
	public HashMap<String, String> m_sampleGroups = new HashMap<String, String>();
	/**a HashMap stores information of all sample groups*/
	public HashMap<String, String> m_groupTypes = new HashMap<String, String>();
	/**an ArrayList of strings stores all fastq files*/
	public ArrayList<String[]> m_fastqList = new  ArrayList<String[]>();
	/**an ArrayList stores all BAM files*/
	public ArrayList<String[]> m_bamList = new  ArrayList<String[]>();
	
	private static YamlParser instance = null;
	
	/**default constructor for the class YamlParse*/
	protected YamlParser() {
		// Exists only to defeat instantiation.
	}
	
	/**
	 * create a new YamlParse instance
	 * @return a YamlParse instance
	 */
	public static YamlParser getInstance() {
		if (instance == null) {
			instance = new YamlParser();
		}
		return instance;
	}
	   
	/**
	 * a function to parse customer's yaml file to provide informations of 
	 * customers' demands 
	 * @param yamlFile a yamlFile that contains the demands of customers
	 * @throws IOException -  I/O operations are failed are interrupted.
	 */
	public void parse(String yamlFile) throws Exception 
	{		
		try
		{	
		    InputStream input = new FileInputStream(new File(yamlFile));
		    Yaml yaml = new Yaml();
		    System.out.println();
		    System.out.println("===========================");
		    System.out.println("Yaml file: " + yamlFile);
            /*print out the yaml file information*/
		    for (Object data : yaml.loadAll(input)) {  //?
                /*a map instance to store the data of the yaml file*/
                Map map = (Map) data;//?
		    		
                /*store the value associated with the key "project" to
                 * m_projectName
                 */
                if (map.get("project") != null)
                {
                    m_projectName = (String) map.get("project");
                    System.out.println("project: " + m_projectName);
                }
                if (map.get("analysis") != null)
                {
                    m_analysis = (String) map.get("analysis");
                    System.out.println("analysis: " + m_analysis);
                }
                if (map.get("date") != null)
                {
                    m_projectDate = (String) map.get("date");
                    System.out.println("date: " + m_projectDate);
                }
                if (map.get("upload") != null)
                {
                    m_S3UploadURL = (String) map.get("upload");
                    System.out.println("upload: " + m_S3UploadURL);
                }
		    		
                /*a map instance to store the data of each sample in the yaml file*/
                Map sample = (Map) map.get("sample");
                if (sample != null)
                {
                    String fileName = (String) sample.get("filename");
                    String S3DownloadURL = (String) sample.get("download");
                    String description = (String) sample.get("description");
                    String group = (String) sample.get("group");
                    
                    System.out.println("---");
                    System.out.println("sample:");
                    System.out.println("  filename: " + fileName);
                    System.out.println("  download: " + S3DownloadURL);
                    System.out.println("  description: " + description);
                    System.out.println("  group: " + group);
			    		
                    //Keep the fastq files information to the m_fastqList
                    if (fileName.indexOf(".fq.gz") > -1 || fileName.indexOf(".fastq") > -1)
                        m_fastqList.add(new String[]{fileName, S3DownloadURL, description, group});

                    //Keep the BAM files information to the m_bamList
                    if (fileName.indexOf(".bam") > -1)
                        m_bamList.add(new String[]{fileName.substring(0, fileName.indexOf(".bam")), S3DownloadURL, description, group});
			    		
					//Keep the group information
					if (group != null && !"NA".equalsIgnoreCase(group) && ! m_sampleGroups.containsKey(fileName))
						m_sampleGroups.put(fileName, group);
					if (group != null && !"NA".equalsIgnoreCase(group) && !m_groupTypes.containsKey(group))
					{
						m_groupTypes.put(group, group);
					}
	    			}
		    }
		    
		    System.out.println("===========================");
		    System.out.println();
		    
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * a function to return the date of the project
	 * @return the data of the project
	 */
	public String getProjectDate()
	{
		return m_projectDate;
	}
	
	/**
	 * a function to return methods of analysis
	 * @return methods of analysis
	 */
	public String getAnalysis()
	{
		return m_analysis;
	}	
	
	/**
	 * a function to return the S3 upload URL
	 * @return S3 upload URL
	 */
	public String getS3UploadURL()
	{
		return m_S3UploadURL;
	}
	
	/**
	 * a function to return all fastq files
	 * @return all fastq files
	 */
	public ArrayList<String[]> getFastqList()
	{
		return m_fastqList;
	}
	
	/**
	 * a function to return all BAM files
	 * @return all BAM files
	 */
	public ArrayList<String[]> getBamList()
	{
		return m_bamList;
	}
	
	/**
	 * a function to return information of all sample groups
	 * @return information of all sample groups
	 */
	public HashMap<String, String> getSampleGroups()
	{
		return m_sampleGroups;
	}
}
