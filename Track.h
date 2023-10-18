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