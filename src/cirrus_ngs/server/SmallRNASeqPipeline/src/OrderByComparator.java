
import java.util.Comparator;

class OrderByComparator implements Comparator<Object> {  
        public int compare(Object o1,Object o2) {  
                Field f1 = (Field)o1;  
                Field f2 = (Field)o2;  
                int i = f1.getFamilyName().compareTo(f2.getFamilyName());  
                int j = f1.getType().compareTo(f2.getType());  
                float k = Float.valueOf(f2.getValue()) - Float.valueOf(f1.getValue());  
                int retStatus=0;  
                if ( i < 0 ) {  
                        retStatus = -1;  
                } else if ( i > 0 ) {  
                        retStatus = 1;  
                } else {  
                        if ( j < 0 ) {  
                                retStatus = -1;  
                        } else if ( j > 0 ) {  
                                retStatus = 1;  
                        } else {  
                                if ( k < 0 ) {  
                                        retStatus = -1;  
                                } else if( k > 0 ){  
                                        retStatus = 1;  
                                } else {  
                                        retStatus = 0;  
                                }  
                        }  
                }  
                return retStatus;  
        }  
}  