
import java.util.*;  
import java.io.*;  

/**
 * The class is to define the all annotation fields
 * @author guorong
 *
 */
class Field {  
        String mirName;  
        String familyName;  
        String isoform;  
        String readID;  
        String type;  
        String value; 
        Field(String mirName,String familyName,String isoform, String readID,String type,String value) {  
                this.mirName = mirName;  
                this.familyName = familyName;  
                this.isoform = isoform;  
                this.readID = readID;  
                this.type = type;  
                this.value = value;  
        }  
        String getMirName() { return mirName; }  
        String getFamilyName() { return familyName; }  
        String getIsoform() { return isoform; }  
        String getReadID() { return readID; }  
        String getType() { return type; }  
        String getValue() { return value; }  
        public String toString() {  
                return mirName+","+familyName+","+isoform+","+readID+","+type+","+value;  
        }  
} 