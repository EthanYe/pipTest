#ifndef HYDROMENT_H
#define HYDROMENT_H

#include <QGraphicsPixmapItem>
#include <QVector>

QT_BEGIN_NAMESPACE
class QPixmap;
class QGraphicsSceneContextMenuEvent;
class QMenu;
class QPolygonF;
QT_END_NAMESPACE

class Qpipe;


//! [0]
class Hydroment : public QGraphicsItem
{
public:

    enum { Type = UserType + 15 };
    enum HydroType { Step, Conditional, StartEnd, Io };

    Hydroment(HydroType hydroType, QMenu *contextMenu, QGraphicsItem *parent = nullptr);

    void removePipe(Pipe *pipe);
    void removePipes();
    HydroType hydroType() const { return myHydroType; }
    QPolygonF polygon() const { return myPolygon; }
    void addPipe(Pipe *Pipe);
    QPixmap image() const;
    int type() const override { return Type; }

protected:
    void contextMenuEvent(QGraphicsSceneContextMenuEvent *event) override;
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override;

private:
    HydroType myHydroType;
    QPolygonF myPolygon;
    QMenu *myContextMenu;
    QVector<Pipe *> Pipes;
};
//! [0]

#endif // HYDROMENT_H
