#ifndef TRACK_H
#define TRACK_H

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

#endif // #define TRACK_H