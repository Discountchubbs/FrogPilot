/****************************************************************************
** Meta object code from reading C++ file 'sidebar.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.8)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "sidebar.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'sidebar.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.8. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Sidebar_t {
    QByteArrayData data[20];
    char stringdata0[200];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Sidebar_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Sidebar_t qt_meta_stringdata_Sidebar = {
    {
QT_MOC_LITERAL(0, 0, 7), // "Sidebar"
QT_MOC_LITERAL(1, 8, 12), // "openSettings"
QT_MOC_LITERAL(2, 21, 0), // ""
QT_MOC_LITERAL(3, 22, 5), // "index"
QT_MOC_LITERAL(4, 28, 5), // "param"
QT_MOC_LITERAL(5, 34, 12), // "valueChanged"
QT_MOC_LITERAL(6, 47, 17), // "offroadTransition"
QT_MOC_LITERAL(7, 65, 7), // "offroad"
QT_MOC_LITERAL(8, 73, 11), // "updateState"
QT_MOC_LITERAL(9, 85, 7), // "UIState"
QT_MOC_LITERAL(10, 93, 1), // "s"
QT_MOC_LITERAL(11, 95, 13), // "connectStatus"
QT_MOC_LITERAL(12, 109, 10), // "ItemStatus"
QT_MOC_LITERAL(13, 120, 11), // "pandaStatus"
QT_MOC_LITERAL(14, 132, 10), // "tempStatus"
QT_MOC_LITERAL(15, 143, 7), // "netType"
QT_MOC_LITERAL(16, 151, 11), // "netStrength"
QT_MOC_LITERAL(17, 163, 9), // "cpuStatus"
QT_MOC_LITERAL(18, 173, 12), // "memoryStatus"
QT_MOC_LITERAL(19, 186, 13) // "storageStatus"

    },
    "Sidebar\0openSettings\0\0index\0param\0"
    "valueChanged\0offroadTransition\0offroad\0"
    "updateState\0UIState\0s\0connectStatus\0"
    "ItemStatus\0pandaStatus\0tempStatus\0"
    "netType\0netStrength\0cpuStatus\0"
    "memoryStatus\0storageStatus"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Sidebar[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       8,   60, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    2,   44,    2, 0x06 /* Public */,
       1,    1,   49,    2, 0x26 /* Public | MethodCloned */,
       1,    0,   52,    2, 0x26 /* Public | MethodCloned */,
       5,    0,   53,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       6,    1,   54,    2, 0x0a /* Public */,
       8,    1,   57,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int, QMetaType::QString,    3,    4,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void, QMetaType::Bool,    7,
    QMetaType::Void, 0x80000000 | 9,   10,

 // properties: name, type, flags
      11, 0x80000000 | 12, 0x0049500b,
      13, 0x80000000 | 12, 0x0049500b,
      14, 0x80000000 | 12, 0x0049500b,
      15, QMetaType::QString, 0x00495003,
      16, QMetaType::Int, 0x00495003,
      17, 0x80000000 | 12, 0x0049500b,
      18, 0x80000000 | 12, 0x0049500b,
      19, 0x80000000 | 12, 0x0049500b,

 // properties: notify_signal_id
       3,
       3,
       3,
       3,
       3,
       3,
       3,
       3,

       0        // eod
};

void Sidebar::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Sidebar *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->openSettings((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< const QString(*)>(_a[2]))); break;
        case 1: _t->openSettings((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->openSettings(); break;
        case 3: _t->valueChanged(); break;
        case 4: _t->offroadTransition((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: _t->updateState((*reinterpret_cast< const UIState(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (Sidebar::*)(int , const QString & );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Sidebar::openSettings)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (Sidebar::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Sidebar::valueChanged)) {
                *result = 3;
                return;
            }
        }
    } else if (_c == QMetaObject::RegisterPropertyMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 7:
        case 6:
        case 5:
        case 2:
        case 1:
        case 0:
            *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< ItemStatus >(); break;
        }
    }

#ifndef QT_NO_PROPERTIES
    else if (_c == QMetaObject::ReadProperty) {
        auto *_t = static_cast<Sidebar *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< ItemStatus*>(_v) = _t->connect_status; break;
        case 1: *reinterpret_cast< ItemStatus*>(_v) = _t->panda_status; break;
        case 2: *reinterpret_cast< ItemStatus*>(_v) = _t->temp_status; break;
        case 3: *reinterpret_cast< QString*>(_v) = _t->net_type; break;
        case 4: *reinterpret_cast< int*>(_v) = _t->net_strength; break;
        case 5: *reinterpret_cast< ItemStatus*>(_v) = _t->cpu_status; break;
        case 6: *reinterpret_cast< ItemStatus*>(_v) = _t->memory_status; break;
        case 7: *reinterpret_cast< ItemStatus*>(_v) = _t->storage_status; break;
        default: break;
        }
    } else if (_c == QMetaObject::WriteProperty) {
        auto *_t = static_cast<Sidebar *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0:
            if (_t->connect_status != *reinterpret_cast< ItemStatus*>(_v)) {
                _t->connect_status = *reinterpret_cast< ItemStatus*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        case 1:
            if (_t->panda_status != *reinterpret_cast< ItemStatus*>(_v)) {
                _t->panda_status = *reinterpret_cast< ItemStatus*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        case 2:
            if (_t->temp_status != *reinterpret_cast< ItemStatus*>(_v)) {
                _t->temp_status = *reinterpret_cast< ItemStatus*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        case 3:
            if (_t->net_type != *reinterpret_cast< QString*>(_v)) {
                _t->net_type = *reinterpret_cast< QString*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        case 4:
            if (_t->net_strength != *reinterpret_cast< int*>(_v)) {
                _t->net_strength = *reinterpret_cast< int*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        case 5:
            if (_t->cpu_status != *reinterpret_cast< ItemStatus*>(_v)) {
                _t->cpu_status = *reinterpret_cast< ItemStatus*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        case 6:
            if (_t->memory_status != *reinterpret_cast< ItemStatus*>(_v)) {
                _t->memory_status = *reinterpret_cast< ItemStatus*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        case 7:
            if (_t->storage_status != *reinterpret_cast< ItemStatus*>(_v)) {
                _t->storage_status = *reinterpret_cast< ItemStatus*>(_v);
                Q_EMIT _t->valueChanged();
            }
            break;
        default: break;
        }
    } else if (_c == QMetaObject::ResetProperty) {
    }
#endif // QT_NO_PROPERTIES
}

QT_INIT_METAOBJECT const QMetaObject Sidebar::staticMetaObject = { {
    &QFrame::staticMetaObject,
    qt_meta_stringdata_Sidebar.data,
    qt_meta_data_Sidebar,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *Sidebar::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Sidebar::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Sidebar.stringdata0))
        return static_cast<void*>(this);
    return QFrame::qt_metacast(_clname);
}

int Sidebar::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QFrame::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
#ifndef QT_NO_PROPERTIES
    else if (_c == QMetaObject::ReadProperty || _c == QMetaObject::WriteProperty
            || _c == QMetaObject::ResetProperty || _c == QMetaObject::RegisterPropertyMetaType) {
        qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 8;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 8;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 8;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 8;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 8;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}

// SIGNAL 0
void Sidebar::openSettings(int _t1, const QString & _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 3
void Sidebar::valueChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 3, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
