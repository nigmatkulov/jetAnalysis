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