/**
 * @file GenTrack.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Generated track description
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef GenTrack_h
#define GenTrack_h

// ROOT headers
#include "TObject.h"

//_________________
class GenTrack : public TObject {
  public:
    GenTrack();
    virtual ~GenTrack();

  private:

    ClassDef(GenTrack, 1)
};

#endif // #define GenTrack_h