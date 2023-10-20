/**
 * @file Track.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Track class description
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef Track_h
#define Track_h

// ROOT headers
#include "TObject.h"

//________________
class Track : public TObject {
  public:
    Track();
    virtual ~Track();

  private:
    
    ClassDef(Track, 1)
};

#endif // #define Track_h