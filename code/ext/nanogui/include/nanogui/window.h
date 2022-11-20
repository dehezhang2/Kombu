#if !defined(__NANOGUI_WINDOW_H)
#define __NANOGUI_WINDOW_H

#include <nanogui/widget.h>

NANOGUI_NAMESPACE_BEGIN

class Window : public Widget {
    friend class Popup;
public:
    Window(Widget *parent, const std::string &title = "Untitled");

    /// Return the window title
    inline const std::string &title() const { return mTitle; }
    /// Set the window title
    inline void setTitle(const std::string &title) { mTitle = title; }

    /// Is this a model dialog?
    inline bool modal() const { return mModal; }
    /// Set whether or not this is a modal dialog
    inline void setModal(bool modal) { mModal = modal; }

    /// Dispose the window
    void dispose();

    /// Center the window in the current \ref Screen
    void center();

    /// Draw the window
    virtual void draw(NVGcontext *ctx);

    /// Handle window drag events
    virtual bool mouseDragEvent(const Vector2i &p, const Vector2i &rel, int button, int modifiers);
    /// Handle mouse events recursively and bring the current window to the top
    virtual bool mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers);
    /// Accept scroll events and propagate them to the widget under the mouse cursor
    virtual bool scrollEvent(const Vector2i &p, const Vector2f &rel);
protected:
    /// Internal helper function to maintain nested window position values; overridden in \ref Popup
    virtual void refreshRelativePlacement();
protected:
    std::string mTitle;
    bool mModal;
};

NANOGUI_NAMESPACE_END

#endif /* __NANOGUI_WINDOW_H */
