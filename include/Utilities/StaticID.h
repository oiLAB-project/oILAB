/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_STATICID_H_
#define gbLAB_STATICID_H_


namespace gbLAB
{
	
	/************************************************************/	
	/************************************************************/
	/*! \brief A class template that implements a counter of the number of 
     * instances of Derived type that are created at runtime. It also provides a 
     * unique increasing static identifier (sID) for each instance.
	 *
     * Example:
     * \include test/test_StaticID/main.cpp 
     * Output:
     * \include test/test_StaticID/output.txt
	 */
	template<typename Derived>
	class StaticID
    {

		// The increment
		static size_t increment;
		
		// The incremental counters
		static size_t count;
        static bool count_used;
		
	public:
		
		//! The static ID of this
		const  size_t sID;
		
        /**********************************************************************/
		StaticID() : sID(count)
        {
            count_used=true;
            count+=increment;
		}
		
        /**********************************************************************/
		StaticID(const StaticID&) : sID(count)
        {
            count_used=true;
            count+=increment;
		}
        
        /**********************************************************************/
        static size_t nextID()
        {
            return count;
        }
        
        static size_t& get_count()
        {
            return count;
        }
		
        /**********************************************************************/
		static void set_count(const size_t& newCount)
        {
            if (newCount<count)
            {
                throw std::runtime_error("StaticID error: YOU ARE TRYING TO SET THE COUNTER TO A LOWER VALUE THAN THE CURRENT ONE\n");
            }
            count =  newCount;
            count_used=false;
		}
		
        /**********************************************************************/
		static void set_increment(const size_t& newIncrement)
        {
            if (newIncrement<1)
            {
                throw std::runtime_error("StaticID error: newIncrement MUST BE >=1\n");
            }
            if(count_used)
            {
                count-=increment;
                count+=newIncrement;
            }
            increment = newIncrement;
		}
		
	};
	
	/* Static data members  *****************************/
	template<typename Derived>
	size_t StaticID<Derived>::increment = 1;

	template<typename Derived>
	size_t StaticID<Derived>::count = 0;

    template<typename Derived>
    bool StaticID<Derived>::count_used = false;
	
} // namespace gbLAB
#endif
