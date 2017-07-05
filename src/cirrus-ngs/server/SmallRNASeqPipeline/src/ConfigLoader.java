import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

/**
 * The class is to load the configuration file which contains all parameters and path of modules
 * the pipeline.
 * @author guorong
 *
 */
public class ConfigLoader {
	public HashMap<String, String>  m_configsList = new HashMap<String, String>();
	
	public void parsePathFile(String file) throws Exception
	{
		String stringLine = null;
		try
		{		
			BufferedReader in = new BufferedReader(new FileReader(file));
			
			while ((stringLine = in.readLine()) != null)
			{
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
	
	private void processPathFile(String stringLine) throws IOException
	{
		StringTokenizer st = new StringTokenizer(stringLine, "\t");
		
		int columnth = 0;
		String tool = null;
		String path = null;
			
		if ((stringLine.length() == 0))
			return;

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
}
