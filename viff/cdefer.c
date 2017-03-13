/*
 * Copyright 2010 VIFF Development Team.
 *
 * This file is part of VIFF, the Virtual Ideal Functionality Framework.
 *
 * VIFF is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License (LGPL) as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * VIFF is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with VIFF. If not, see <http://www.gnu.org/licenses/>.
 *
 * Original C implementation of Deferred by
 * http://twistedmatrix.com/trac/ticket/2245
 *
 */

#include <Python.h>
#include "structmember.h"

/* Py_VISIT and Py_CLEAR are defined here to be compatible with Python 2.3 */

#ifndef Py_VISIT
#define Py_VISIT(op) \
    do { \
        if (op) { \
            int vret = visit((PyObject *)(op), arg); \
            if (vret) \
                return vret; \
        } \
    } while (0)
#endif

#ifndef Py_CLEAR
#define Py_CLEAR(op) \
    do { \
        if (op) { \
            PyObject *tmp = (PyObject *)(op); \
            (op) = NULL; \
            Py_DECREF(tmp); \
        } \
    } while (0)
#endif

/* be compatible with Python < 2.6 */

#ifndef Py_TYPE
#define Py_TYPE(ob)             (((PyObject*)(ob))->ob_type)
#endif

PyObject * failure_class = NULL;
PyObject * already_called = NULL;
PyObject * debuginfo_class = NULL;
PyObject * format_stack = NULL;

struct cdefer_Callback_ {
    PyObject *func;
    PyObject *args;
    PyObject *kwargs;
    struct cdefer_Callback_ *next;
};

typedef struct cdefer_Callback_ cdefer_Callback;

typedef struct {
    PyObject_HEAD
    PyObject *result;
    int paused;
    cdefer_Callback *callbacks;
    cdefer_Callback *last_callback;
    PyObject *debuginfo;
    int called;
    int runningCallbacks;
} cdefer_Deferred;

/* Prototypes */

static PyObject *cdefer_setDebugging(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs);

static PyObject *cdefer_getDebugging(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs);

static PyObject * cdefer_Deferred_new(PyTypeObject *type, PyObject *args,
        PyObject *kwargs);

static void cdefer_Deferred_dealloc(PyObject *o);

static int cdefer_Deferred_traverse(PyObject *o, visitproc visit, void *arg);

static int cdefer_Deferred_clear(PyObject *o);

static int cdefer_Deferred_clear(PyObject *o);

static int cdefer_Deferred___init__(cdefer_Deferred *self, PyObject *args,
        PyObject *kwargs);

static PyObject *cdefer_Deferred__addCallbacks(cdefer_Deferred *self,
        PyObject *callback, PyObject *errback, PyObject *callbackArgs,
        PyObject *callbackKeywords, PyObject *errbackArgs,
        PyObject *errbackKeywords);

static PyObject *cdefer_Deferred_addCallback(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs);

static PyObject *cdefer_Deferred_addErrback(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs);

static PyObject *cdefer_Deferred_addBoth(cdefer_Deferred *self, PyObject *args,
        PyObject *kwargs);

static PyObject *cdefer_Deferred_pause(cdefer_Deferred *self, PyObject *args);

static PyObject *cdefer_Deferred_unpause(cdefer_Deferred *self,
        PyObject *args);

static PyObject *cdefer_Deferred_chainDeferred(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs);

static PyObject *cdefer_Deferred__runCallbacks(cdefer_Deferred *self);

static PyObject *cdefer_Deferred__startRunCallbacks(cdefer_Deferred *self,
        PyObject *result);

static PyObject *cdefer_Deferred_callback(cdefer_Deferred *self, PyObject *args,
        PyObject *kwargs);

static PyObject *cdefer_Deferred_errback(cdefer_Deferred *self, PyObject *args,
        PyObject *kwargs);

static PyObject *cdefer_Deferred__continue(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs);


static int is_debug = 0;

static PyObject *cdefer_setDebugging(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs)
{
    int new_debug;
    PyObject *on;
    static char *argnames[] = {"on", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", argnames, &on)) {
        return NULL;
    }
    new_debug = PyObject_IsTrue(on);
    if (-1 == new_debug) {
        return NULL;
    }
    is_debug = new_debug;
    Py_INCREF(Py_None);
    return Py_None;
}

static char cdefer_setDebugging_doc[] = "Enable or disable Deferred debugging.\n\n    When debugging is on, the call stacks from creation and invocation are\n    recorded, and added to any AlreadyCalledErrors we raise.\n";


static PyObject *cdefer_getDebugging(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs)
{
    static char *argnames[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "", argnames)) {
        return NULL;
    }
    return PyBool_FromLong(is_debug);
}

static char cdefer_getDebugging_doc[] = "Determine whether Deferred debugging is enabled.\n";


static PyTypeObject cdefer_DeferredType;

static PyObject * cdefer_Deferred_new(PyTypeObject *type, PyObject *args,
                                      PyObject *kwargs)
{
    cdefer_Deferred *self;
    self = (cdefer_Deferred *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static cdefer_Callback* cdefer_Callback_free(cdefer_Callback *callback) {
    cdefer_Callback* next = callback->next;
    Py_DECREF(callback->func);
    Py_DECREF(callback->args);
    Py_DECREF(callback->kwargs);
    free(callback);
    return next;
}

static void cdefer_Deferred_dealloc(PyObject *o) {
    cdefer_Deferred *self;
    self = (cdefer_Deferred *)o;
    PyObject_GC_UnTrack(self);
    Py_XDECREF(self->result);
    Py_XDECREF(self->debuginfo);
    while (self->callbacks)
	self->callbacks = cdefer_Callback_free(self->callbacks);
    (*o->ob_type->tp_free)(o);
}

static int cdefer_Deferred_traverse(PyObject *o, visitproc visit, void *arg) {
    cdefer_Deferred *self;
    cdefer_Callback *cb;
    self = (cdefer_Deferred *)o;
    Py_VISIT(self->result);
    Py_VISIT(self->debuginfo);
    cb = self->callbacks;
    while (cb) {
	Py_VISIT(cb->func);
	Py_VISIT(cb->args);
	Py_VISIT(cb->kwargs);
	cb = cb->next;
    }
    return 0;
}

static int cdefer_Deferred_clear(PyObject *o) {
    cdefer_Deferred *self;
    cdefer_Callback *next;
    self = (cdefer_Deferred *)o;
    Py_CLEAR(self->result);
    Py_CLEAR(self->debuginfo);
    while (self->callbacks) {
	Py_CLEAR(self->callbacks->func);
	Py_CLEAR(self->callbacks->args);
	Py_CLEAR(self->callbacks->kwargs);
	next = self->callbacks->next;
	free(self->callbacks);
	self->callbacks = next;
    }
    return 0;
}

static int cdefer_Deferred__set_debug_stack(cdefer_Deferred *self, char *name)
{
    int rc;
    PyObject *stack;

    /* Keep the debug info object even if we fail to format stack
     * or place it into the dict. */
    stack = PyObject_CallObject(format_stack, NULL);
    if (!stack) {
        return -1;
    }
    rc = PyObject_SetAttrString(self->debuginfo, name, stack);
    /* Unlike other functions of this naming convention (and
     * unlike PyDict_GetItemString), PyDict_SetItemString
     * copies/creates a new reference, so we shouldn't keep ours
     * too. */
    Py_DECREF(stack);
    if (-1 == rc) {
        return -1;
    }
    return 0;
}

static int cdefer_Deferred___init__(cdefer_Deferred *self, PyObject *args,
                                    PyObject *kwargs)
{
    /*
    static char *argnames[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "", argnames)) {
        return -1;
    }
    */
    if (is_debug) {
        self->debuginfo = PyObject_CallObject(debuginfo_class, NULL);
        if (!self->debuginfo) {
            return -1;
        }
        if (-1 == cdefer_Deferred__set_debug_stack(self, "creator")) {
            return -1;
        }
    }

    self->paused = 0;
    self->called = 0;
    self->runningCallbacks = 0;
    self->callbacks = NULL;
    self->last_callback = NULL;
    return 0;
}

static PyObject* cdefer_Deferred___str__(cdefer_Deferred *self) {
    PyObject *result, *tmp;
    const char* cname = Py_TYPE(self)->tp_name;

    if (self->result) {
	tmp = PyObject_Repr(self->result);

	if (!tmp)
	    return NULL;

	result = PyString_FromFormat("<%s at %p  current result: %s>",
				     cname, self, PyString_AsString(tmp));
	Py_DECREF(tmp);
	return result;
    }

    return PyString_FromFormat("<%s at %p>", cname, self);
}

static PyObject *cdefer_Deferred__addCallbacks(cdefer_Deferred *self,
        PyObject *callback, PyObject *errback, PyObject *callbackArgs,
        PyObject *callbackKeywords, PyObject *errbackArgs,
        PyObject *errbackKeywords) {
    PyObject *result;
    cdefer_Callback *cb;

    if (!PyCallable_Check(callback)) {
	PyErr_SetNone(PyExc_AssertionError);
	return NULL;
    }

    /*
    if (callback != Py_None) {
        if (!PyCallable_Check(callback)) {
            PyErr_SetNone(PyExc_AssertionError);
            return NULL;
        }
    }
    if (errback != Py_None) {
        if (!PyCallable_Check(errback)) {
            PyErr_SetNone(PyExc_AssertionError);
            return NULL;
        }
    }
    */

    cb = (cdefer_Callback*) malloc(sizeof(cdefer_Callback));
    cb->next = NULL;
    cb->func = callback;
    cb->args = callbackArgs;
    cb->kwargs = callbackKeywords;
    Py_INCREF(callback);
    Py_INCREF(callbackArgs);
    Py_INCREF(callbackKeywords);

    if (self->last_callback)
	self->last_callback->next = cb;
    else
	self->callbacks = cb;

    self->last_callback = cb;

    if (self->called) {
        if (cdefer_Deferred__runCallbacks(self) == NULL) {
            return NULL;
        }
    }

    result = (PyObject *)self;
    Py_INCREF(result);
    return result;
}

static char cdefer_Deferred_addCallbacks_doc[] = "Add a pair of callbacks (success and error) to this Deferred.\n\nThese will be executed when the \'master\' callback is run.";

static PyObject *cdefer_Deferred_addCallbacks(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs) {
    static char *argnames[] = {"callback", "errback", "callbackArgs",
        "callbackKeywords", "errbackArgs", "errbackKeywords", NULL};
    PyObject *callback;
    PyObject *errback = Py_None;
    PyObject *callbackArgs = Py_None;
    PyObject *callbackKeywords = Py_None;
    PyObject *errbackArgs = Py_None;
    PyObject *errbackKeywords = Py_None;
    PyObject *result;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|OOOOO", argnames,
                &callback, &errback, &callbackArgs,
                &callbackKeywords, &errbackArgs, &errbackKeywords)) {
        return NULL;
    }
    result = cdefer_Deferred__addCallbacks(self, callback, errback,
        callbackArgs, callbackKeywords, errbackArgs, errbackKeywords);
    return result;
}

static char cdefer_Deferred_addCallback_doc[] = "Convenience method for adding just a callback.\n\nSee L{addCallbacks}.";

/* Returns a NEW reference to the callback/errback arg, and to the
 * cbackArgs, but a BORROWED reference to the keywords. In case of
 * error, no references are returned/touched */
static PyObject *extract_cback_args_kw(char *argname,
                                       PyObject *args, PyObject *kwargs,
                                       PyObject **cbackArgs,
                                       PyObject **cbackKeywords)
{
    PyObject *cback;

    if (kwargs) {
        (*cbackKeywords) = kwargs;
    } else {
        (*cbackKeywords) = Py_None;
    }
    if (PyTuple_Size(args) > 0) {
        cback = PyTuple_GET_ITEM(args, 0);
        if (!cback) {
            return NULL;
        }
        (*cbackArgs) = PyTuple_GetSlice(args, 1, PyTuple_Size(args));
        if (!(*cbackArgs)) {
            return NULL;
        }
        Py_INCREF(cback);
    } else {
        cback = PyDict_GetItemString((*cbackKeywords), argname);
        if (!cback) {
            PyErr_Format(PyExc_TypeError, "addCallback requires '%s' argument'", argname);
            return NULL;
        }
        (*cbackArgs) = Py_None;
        Py_INCREF(Py_None);

        /* "callback" in the keyword dict may be the only reference to
         * it, and we delete it from the dict, so we must own a
         * reference too */
        Py_INCREF(cback);

        if (PyDict_DelItemString((*cbackKeywords), argname) == -1) {
            Py_DECREF(cback);
            Py_DECREF(Py_None);
            return NULL;
        }
    }
    return cback;
}

static PyObject *cdefer_Deferred_addCallback(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs) {
    PyObject *callback;
    PyObject *callbackArgs;
    PyObject *callbackKeywords;
    PyObject *result;
    callback = extract_cback_args_kw("callback", args, kwargs, &callbackArgs, &callbackKeywords);
    if (!callback) {
        return NULL;
    }
    result = cdefer_Deferred__addCallbacks(self, callback, Py_None, callbackArgs,
        callbackKeywords, Py_None, Py_None);
    Py_DECREF(callback);
    Py_DECREF(callbackArgs);
    return result;
}

static char cdefer_Deferred_addErrback_doc[] = "Convenience method for adding just an errback.\n\nSee L{addCallbacks}.";

static PyObject *cdefer_Deferred_addErrback(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs) {
    PyObject *errback;
    PyObject *errbackArgs;
    PyObject *errbackKeywords;
    PyObject *result;
    errback = extract_cback_args_kw("errback", args, kwargs, &errbackArgs, &errbackKeywords);
    if (!errback) {
        return NULL;
    }
    /*
    result = cdefer_Deferred__addCallbacks(self, Py_None, errback, Py_None,
        Py_None, errbackArgs, errbackKeywords);
    */
    Py_DECREF(errback);
    Py_DECREF(errbackArgs);

    Py_INCREF(self);
    result = (PyObject *)self;
    return result;
}

static char cdefer_Deferred_addBoth_doc[] = "Convenience method for adding a single callable as both a callback\nand an errback.\n\nSee L{addCallbacks}.";

static PyObject *cdefer_Deferred_addBoth(cdefer_Deferred *self, PyObject *args,
        PyObject *kwargs) {
    PyObject *callback;
    PyObject *callbackArgs;
    PyObject *callbackKeywords;
    PyObject *result;
    callback = extract_cback_args_kw("callback", args, kwargs, &callbackArgs, &callbackKeywords);
    if (!callback) {
        return NULL;
    }
    result = cdefer_Deferred__addCallbacks(self, callback, callback,
        callbackArgs, callbackKeywords, callbackArgs, callbackKeywords);
    Py_DECREF(callback);
    Py_DECREF(callbackArgs);
    return result;
}

static char cdefer_Deferred_pause_doc[] = "Stop processing on a Deferred until L{unpause}() is called.";

static PyObject *cdefer_Deferred_pause(cdefer_Deferred *self, PyObject *args) {
    PyObject *result;
    self->paused++;
    result = Py_None;
    Py_INCREF(Py_None);
    return result;
}

static char cdefer_Deferred_unpause_doc[] = "Process all callbacks made since L{pause}() was called.";

static PyObject *cdefer_Deferred_unpause(cdefer_Deferred *self,
        PyObject *args) {
    self->paused--;
    if (!self->paused && self->called) {
        return cdefer_Deferred__runCallbacks(self);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static char cdefer_Deferred_chainDeferred_doc[] = "Chain another Deferred to this Deferred.\n\nThis method adds callbacks to this Deferred to call d\'s callback or\nerrback, as appropriate.";

static PyObject *cdefer_Deferred_chainDeferred(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs) {
    PyObject *d;
    PyObject *callback;
    PyObject *errback;
    PyObject *result;
    static char *argnames[] = {"d", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", argnames, &d)) {
        return NULL;
    }
    callback = PyObject_GetAttrString(d, "callback");
    if (!callback) {
        return NULL;
    }
    errback = PyObject_GetAttrString(d, "errback");
    if (!errback) {
        Py_DECREF(callback);
        return NULL;
    }
    result = cdefer_Deferred__addCallbacks(self, callback, errback, Py_None,
        Py_None, Py_None, Py_None);
    Py_DECREF(callback);
    Py_DECREF(errback);
    return result;
}

/*
static int cdefer_Deferred__set_debuginfo_fail_result(cdefer_Deferred *self) {
    if (!self->debuginfo) {
        self->debuginfo = PyObject_CallObject(debuginfo_class, NULL);
        if (!self->debuginfo) {
            return -1;
        }
    }
    if (PyObject_SetAttrString(self->debuginfo, "failResult", self->result) == -1) {
        return -1;
    }
    return 0;
}

static int cdefer_Deferred__clear_debuginfo(cdefer_Deferred *self) {
    if (self->debuginfo) {
        if (PyObject_SetAttrString(self->debuginfo, "failResult", Py_None) == -1) {
            return -1;
        }
    }
    return 0;
}
*/

static PyObject *cdefer_Deferred__runCallbacks(cdefer_Deferred *self) {
    PyObject *callback;
    PyObject *args;
    PyObject *newArgs;
    PyObject *newArgs2;
    PyObject *kwargs;
    PyObject *_continue;
    //PyObject *type, *value, *traceback;
    PyObject *tmp;
    PyObject *result;

    if (self->runningCallbacks) {
	Py_INCREF(Py_None);
	return Py_None;
    }

    if (!self->paused) {
        while(self->callbacks) {
	    /*
            if (PyObject_IsInstance(self->result, failure_class)) {
            }
	    */

	    callback = self->callbacks->func;
            if (!callback) {
		PyErr_SetString(PyExc_RuntimeError, "NULL callback.");
                return NULL;
            }
	    /*
            if (callback == Py_None) {
                self->callbacks = cdefer_Callback_free(self->callbacks);
		if (!self->callbacks)
		    self->last_callback = NULL;
                continue;
            }
	    */

	    args = self->callbacks->args;
            if (!args) {
		PyErr_SetString(PyExc_RuntimeError, "NULL args.");
                return NULL;
            }

	    kwargs = self->callbacks->kwargs;
            if (!kwargs) {
		PyErr_SetString(PyExc_RuntimeError, "NULL kwargs.");
                return NULL;
            }

            newArgs = Py_BuildValue("(O)", self->result);
            if (!newArgs) {
                return NULL;
            }

            if (args != Py_None) {
                newArgs2 = PySequence_InPlaceConcat(newArgs, args);
                if (!newArgs2) {
                    return NULL;
                }
                Py_CLEAR(newArgs);
            } else {
                newArgs2 = newArgs;
                newArgs = NULL;
            }

	    self->runningCallbacks = 1;

            if (kwargs == Py_None) {
                tmp = PyObject_Call(callback, newArgs2, NULL);
            } else {
                tmp = PyObject_Call(callback, newArgs2, kwargs);
            }

	    self->runningCallbacks = 0;

            Py_DECREF(self->result);
            self->result = tmp;

            Py_CLEAR(newArgs2);

	    self->callbacks = cdefer_Callback_free(self->callbacks);

	    if (!self->callbacks)
		self->last_callback = NULL;

            if (!self->result) {
		return NULL;
		/*
                PyErr_Fetch(&type, &value, &traceback);
                PyErr_NormalizeException(&type, &value, &traceback);
                if (!traceback) {
                    traceback = Py_None;
                    Py_INCREF(traceback);
                }

                self->result = PyObject_CallFunction(failure_class, "OOO", value, type, traceback);
                if (!self->result) {
                    PyErr_Restore(type, value, traceback);
                    return NULL;
                }
                Py_DECREF(type);
                Py_DECREF(value);
                Py_DECREF(traceback);
                continue;
		*/
            }

            if (PyObject_TypeCheck(self->result, &cdefer_DeferredType)) {
                result = cdefer_Deferred_pause(self, NULL);
                if (!result) {
                    return NULL;
                }
                Py_DECREF(result);

                _continue = PyObject_GetAttrString((PyObject *)self,
                                                   "_continue");
                if (!_continue) {
                    return NULL;
                }

		// self->result might be replaced and thus decref'd
		// while running its callback.
		tmp = self->result;
		Py_INCREF(tmp);
                result = cdefer_Deferred__addCallbacks(
                    (cdefer_Deferred *)self->result, _continue,
                    _continue, Py_None, Py_None, Py_None, Py_None);
                /* The reference was either copied/incref'd or not
                 * (when errored) in addCallbacks, either way, we own
                 * one too, and don't need it anymore. */
                Py_DECREF(_continue);
		Py_DECREF(tmp);

                if (!result) {
                    return NULL;
                }
                Py_DECREF(result);

                goto endLabel;
            }
        }
    }
endLabel:;
    /*
    if (PyObject_IsInstance(self->result, failure_class)) {
        result = PyObject_CallMethod((PyObject *)self->result,
                                     "cleanFailure", NULL);
        if (!result) {
            return NULL;
        }
        Py_DECREF(result);
        if (cdefer_Deferred__set_debuginfo_fail_result(self) == -1) {
            return NULL;
        }
    } else {
        if (cdefer_Deferred__clear_debuginfo(self) == -1) {
            return NULL;
        }
    }
    */
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *cdefer_Deferred__startRunCallbacks(cdefer_Deferred *self,
                                                    PyObject *result)
{
    PyObject * already_called_instance;
    PyObject * debug_tracebacks;

    if (is_debug && !self->debuginfo) {
        self->debuginfo = PyObject_CallObject(debuginfo_class, NULL);
        if (!self->debuginfo) {
            return NULL;
        }
    }

    if (self->called) {
        if (is_debug) {
            debug_tracebacks = PyObject_CallMethod(
                self->debuginfo, "_getDebugTracebacks", "s", "\n");
            if (!debug_tracebacks) {
                return NULL;
            }
            already_called_instance = PyObject_CallFunction(already_called, "O", debug_tracebacks);
            Py_DECREF(debug_tracebacks);
            if (!already_called_instance) {
                return NULL;
            }
            PyErr_SetObject(already_called, already_called_instance);
            Py_DECREF(already_called_instance);
            return NULL;
        }
        PyErr_SetNone(already_called);
        return NULL;
    }
    if (is_debug) {
        if (-1 == cdefer_Deferred__set_debug_stack(self, "invoker")) {
            return NULL;
        }
    }

    self->called = 1;
    Py_XDECREF(self->result);
    self->result = result;
    Py_INCREF(self->result);
    return cdefer_Deferred__runCallbacks(self);
}

static char cdefer_Deferred_callback_doc[] = "Run all success callbacks that have been added to this Deferred.\n\nEach callback will have its result passed as the first\nargument to the next; this way, the callbacks act as a\n\'processing chain\'. Also, if the success-callback returns a Failure\nor raises an Exception, processing will continue on the *error*-\ncallback chain.";

static PyObject *cdefer_Deferred_callback(cdefer_Deferred *self, PyObject *args,
        PyObject *kwargs) {
    PyObject *result;
    static char *argnames[] = {"result", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", argnames, &result)) {
        return NULL;
    }
    return cdefer_Deferred__startRunCallbacks(self, result);
}

static char cdefer_Deferred_errback_doc[] = "Run all error callbacks that have been added to this Deferred.\n\nEach callback will have its result passed as the first\nargument to the next; this way, the callbacks act as a\n\'processing chain\'. Also, if the error-callback returns a non-Failure\nor doesn\'t raise an Exception, processing will continue on the\n*success*-callback chain.\n\nIf the argument that\'s passed to me is not a Failure instance,\nit will be embedded in one. If no argument is passed, a Failure\ninstance will be created based on the current traceback stack.\n\nPassing a string as `fail\' is deprecated, and will be punished with\na warning message.";

static PyObject *cdefer_Deferred_errback(cdefer_Deferred *self, PyObject *args,
        PyObject *kwargs) {
    PyObject *fail;
    PyObject *tmp;
    PyObject *result;
    static char *argnames[] = {"fail", NULL};
    fail = Py_None;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|O", argnames, &fail)) {
        return NULL;
    }

    if (PyObject_IsInstance(fail, failure_class)) {
        /* Make "fail" belong to us even if we don't create a Failure
         * wrapper (If we do, the wrapper belongs to us) */
        Py_INCREF(fail);
    } else {
        tmp = PyObject_CallFunction(failure_class, "(O)", fail);
        if (!tmp) {
            return NULL;
        }
        fail = tmp;
    }
    result = cdefer_Deferred__startRunCallbacks(self, fail);
    Py_DECREF(fail);
    return result;
}

static PyObject *cdefer_Deferred__continue(cdefer_Deferred *self,
        PyObject *args, PyObject *kwargs) {
    PyObject *result;
    static char *argnames[] = {"result", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", argnames, &result)) {
        return NULL;
    }
    Py_XDECREF(self->result);
    self->result = result;
    Py_INCREF(self->result);
    return cdefer_Deferred_unpause(self, NULL);
}

static struct PyMethodDef cdefer_Deferred_methods[] = {
  {"addCallbacks", (PyCFunction)cdefer_Deferred_addCallbacks,
                   METH_VARARGS|METH_KEYWORDS, cdefer_Deferred_addCallbacks_doc},
  {"addCallback", (PyCFunction)cdefer_Deferred_addCallback,
                  METH_VARARGS|METH_KEYWORDS, cdefer_Deferred_addCallback_doc},
  {"addErrback", (PyCFunction)cdefer_Deferred_addErrback,
                 METH_VARARGS|METH_KEYWORDS, cdefer_Deferred_addErrback_doc},
  {"addBoth", (PyCFunction)cdefer_Deferred_addBoth,
               METH_VARARGS|METH_KEYWORDS, cdefer_Deferred_addBoth_doc},
  {"chainDeferred", (PyCFunction)cdefer_Deferred_chainDeferred,
                    METH_VARARGS|METH_KEYWORDS, cdefer_Deferred_chainDeferred_doc},
  {"callback", (PyCFunction)cdefer_Deferred_callback,
               METH_VARARGS|METH_KEYWORDS, cdefer_Deferred_callback_doc},
  {"errback", (PyCFunction)cdefer_Deferred_errback,
              METH_VARARGS|METH_KEYWORDS, cdefer_Deferred_errback_doc},
  {"pause", (PyCFunction)cdefer_Deferred_pause,
            METH_VARARGS, cdefer_Deferred_pause_doc},
  {"unpause", (PyCFunction)cdefer_Deferred_unpause,
              METH_VARARGS, cdefer_Deferred_unpause_doc},
  {"_continue", (PyCFunction)cdefer_Deferred__continue,
                METH_VARARGS|METH_KEYWORDS, ""},
  {0, 0, 0, 0}
};

static struct PyMemberDef cdefer_Deferred_members[] = {
  {"result", T_OBJECT_EX, offsetof(cdefer_Deferred, result), 0, 0},
  {"paused", T_INT, offsetof(cdefer_Deferred, paused), READONLY, 0},
  {"called", T_INT, offsetof(cdefer_Deferred, called), READONLY, 0},
  {0, 0, 0, 0, 0}
};

static PyTypeObject cdefer_DeferredType = {
    PyObject_HEAD_INIT(0)
    0,                          /*ob_size*/
    "cdefer.Deferred",          /*tp_name*/
    sizeof(cdefer_Deferred),    /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)cdefer_Deferred_dealloc,    /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    (reprfunc)cdefer_Deferred___str__,    /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    (reprfunc)cdefer_Deferred___str__,    /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE|Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    "This is a callback which will be put off until later.\n\nWhy do we want this? Well, in cases where a function in a threaded\nprogram would block until it gets a result, for Twisted it should\nnot block. Instead, it should return a Deferred.\n\nThis can be implemented for protocols that run over the network by\nwriting an asynchronous protocol for twisted.internet. For methods\nthat come from outside packages that are not under our control, we use\nthreads (see for example L{twisted.enterprise.adbapi}).\n\nFor more information about Deferreds, see doc/howto/defer.html or\nU{http://www.twistedmatrix.com/documents/howto/defer}.", /*tp_doc*/
    (traverseproc)cdefer_Deferred_traverse,   /*tp_traverse*/
    (inquiry)cdefer_Deferred_clear,           /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    cdefer_Deferred_methods,    /*tp_methods*/
    cdefer_Deferred_members,    /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    (initproc)cdefer_Deferred___init__,   /*tp_init*/
    0,                          /*tp_alloc*/
    cdefer_Deferred_new,        /*tp_new*/
    PyObject_GC_Del,            /*tp_free*/
    0,                          /*tp_is_gc*/
    0,                          /*tp_bases*/
    0,                          /*tp_mro*/
    0,                          /*tp_cache*/
    0,                          /*tp_subclasses*/
    0,                          /*tp_weaklist*/
};

static PyMethodDef cdefer_methods[] = {
    {"setDebugging", (PyCFunction)cdefer_setDebugging,
     METH_VARARGS|METH_KEYWORDS, cdefer_setDebugging_doc},

    {"getDebugging", (PyCFunction)cdefer_getDebugging,
     METH_VARARGS|METH_KEYWORDS, cdefer_getDebugging_doc},

    {NULL}  /* Sentinel */
};

/* VIFF */

typedef struct {
    cdefer_Deferred deferred;
    PyObject *runtime;
    PyObject *field;
	long fp;
} cdefer_Share;

static struct PyMemberDef cdefer_Share_members[] = {
  {"runtime", T_OBJECT_EX, offsetof(cdefer_Share, runtime), 0, 0},
  {"field", T_OBJECT_EX, offsetof(cdefer_Share, field), 0, 0},
//  {"fp", T_INT, offsetof(cdefer_Share, fp), 0, 0},
  {0, 0, 0, 0, 0}
};

static int cdefer_Share__init(cdefer_Share *self, PyObject* runtime,
			      PyObject* field, PyObject* value, long fp);

static int cdefer_Share_init(cdefer_Share *self, PyObject *args,
			     PyObject *kwargs)
{
    long tfp = 0;
    PyObject *runtime, *field, *value = NULL, *fp = NULL;
    static char *argnames[] = {"runtime", "field", "value", "fp", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|O|O", argnames,
				     &runtime, &field, &value, &fp))
        return -1;
    
    if (field == Py_None) {
	PyErr_SetString(PyExc_AssertionError,
			"Cannot construct share without a field.");
	return -1;
    }

    if (!PyCallable_Check(field)) {
	PyErr_SetString(PyExc_AssertionError,
			"The field is not callable, wrong argument?");
	return -1;
    }

	if(fp != NULL){
		if(PyInt_Check(fp)) {
		    if(PyInt_AS_LONG(fp) > 0){
			   tfp = (long)1;
		    }
		}
    }		
    if(tfp != (long)1){
		tfp = (long)0;
	}
    return cdefer_Share__init(self, runtime, field, value, tfp);
}

static int cdefer_Share__init(cdefer_Share *self, PyObject* runtime,
			      PyObject* field, PyObject* value, long fp)
{
    PyObject *result, *tmp;

    Py_INCREF(Py_None);
    if (cdefer_Deferred___init__(&self->deferred, Py_None, Py_None) < 0)
	return -1;
    Py_DECREF(Py_None);

    tmp = self->runtime;
    Py_INCREF(runtime);
    self->runtime = runtime;
    Py_XDECREF(tmp);

    tmp = self->field;
    Py_INCREF(field);
    self->field = field;
    Py_XDECREF(tmp);
    
    self->fp = fp;
    if ((value != Py_None) && (value!= NULL)) {
	result = cdefer_Deferred__startRunCallbacks(&self->deferred, value);
	if (result)
	    Py_DECREF(result);
	else
	    return -1;
    }
    return 0;
}

static void cdefer_Share_dealloc(PyObject *o) {
    cdefer_Share *self;
    self = (cdefer_Share *)o;
    Py_XDECREF(self->runtime);
    Py_XDECREF(self->field);
    cdefer_Deferred_dealloc((PyObject *)&self->deferred);
}

static int cdefer_Share_traverse(PyObject *o, visitproc visit, void *arg) {
    cdefer_Share *self;
    self = (cdefer_Share *)o;
    Py_VISIT(self->runtime);
    Py_VISIT(self->field);
    return cdefer_Deferred_traverse((PyObject *)&self->deferred, visit, arg);
}

static int cdefer_Share_clear(PyObject *o) {
    cdefer_Share *self;
    self = (cdefer_Share *)o;
    Py_CLEAR(self->runtime);
    Py_CLEAR(self->field);
    return cdefer_Deferred_clear((PyObject *)&self->deferred);
}

static PyTypeObject cdefer_ShareType;

static PyObject *cdefer_split_result(PyObject *noarg,
				     PyObject *args) {
    PyObject *result, *clone, *res;

    if(PyArg_UnpackTuple(args, "split_result", 2, 2, &result, &clone) < 0)
	return NULL;
	
	res = cdefer_Deferred__startRunCallbacks((cdefer_Deferred*)clone, result);

    if (!res)
	return NULL;

    Py_DECREF(res);
    Py_INCREF(result);
    return result;
}

PyMethodDef split_result_def = {"split_result",
				(PyCFunction)cdefer_split_result,
				METH_VARARGS, 0};

static PyObject *cdefer_Share_clone(cdefer_Share *self, PyObject *noarg) {
    PyObject *clone, *result, *split_result, *arg;

    Py_INCREF(Py_None);
    clone = PyObject_CallFunctionObjArgs((PyObject *)&cdefer_ShareType,
					 self->runtime, self->field, Py_None, PyInt_FromLong(self->fp), NULL);
    Py_DECREF(Py_None);

    if (!clone)
	return NULL;

    split_result = PyCFunction_New(&split_result_def, NULL);

    if (!split_result)
	return NULL;

    arg = Py_BuildValue("(O)", clone);

    if (!arg)
	return NULL;

    Py_INCREF(Py_None);
    result = cdefer_Deferred__addCallbacks(&(self->deferred),
					   split_result, Py_None,
					   arg, Py_None,
					   Py_None, Py_None);
    Py_DECREF(Py_None);
    Py_DECREF(arg);
    Py_DECREF(split_result);

    if (!result)
	return NULL;

    Py_DECREF(result);
    return clone;
}

static PyObject *cdefer_Share_clone_nofp(cdefer_Share *self, PyObject *noarg) {
    PyObject *clone, *result, *split_result, *arg;

    Py_INCREF(Py_None);
    clone = PyObject_CallFunctionObjArgs((PyObject *)&cdefer_ShareType,
					 self->runtime, self->field, Py_None, PyInt_FromLong(0), NULL);
    Py_DECREF(Py_None);

    if (!clone)
	return NULL;

    split_result = PyCFunction_New(&split_result_def, NULL);

    if (!split_result)
	return NULL;

    arg = Py_BuildValue("(O)", clone);

    if (!arg)
	return NULL;

    Py_INCREF(Py_None);
    result = cdefer_Deferred__addCallbacks(&(self->deferred),
					   split_result, Py_None,
					   arg, Py_None,
					   Py_None, Py_None);
    Py_DECREF(Py_None);
    Py_DECREF(arg);
    Py_DECREF(split_result);

    if (!result)
	return NULL;

    Py_DECREF(result);
    return clone;
}

static PyObject *cdefer_Share_clone_fp(cdefer_Share *self, PyObject *noarg) {
    PyObject *clone, *result, *split_result, *arg;

    Py_INCREF(Py_None);
    clone = PyObject_CallFunctionObjArgs((PyObject *)&cdefer_ShareType,
					 self->runtime, self->field, Py_None, PyInt_FromLong(1), NULL);
    Py_DECREF(Py_None);
    if (!clone)
	return NULL;

    split_result = PyCFunction_New(&split_result_def, NULL);

    if (!split_result)
	return NULL;

    arg = Py_BuildValue("(O)", clone);

    if (!arg)
	return NULL;

    Py_INCREF(Py_None);
    result = cdefer_Deferred__addCallbacks(&(self->deferred),
					   split_result, Py_None,
					   arg, Py_None,
					   Py_None, Py_None);
    Py_DECREF(Py_None);
    Py_DECREF(arg);
    Py_DECREF(split_result);

    if (!result)
	return NULL;

    Py_DECREF(result);
    return clone;
}

static struct PyMethodDef cdefer_Share_methods[] = {
    {"clone", (PyCFunction)cdefer_Share_clone, METH_NOARGS, 0},
    {"clone_nofp", (PyCFunction)cdefer_Share_clone_nofp, METH_NOARGS, 0},
    {"clone_fp", (PyCFunction)cdefer_Share_clone_fp, METH_NOARGS, 0},
    {0, 0, 0, 0}
};

#define SYM_BIN_OP(FUNCNAME, NBNAME, SHORTNAME)	\
static PyObject *FUNCNAME(PyObject *self, PyObject *other) { \
    PyObject* result; \
    PyObject* name = PyString_FromString(#SHORTNAME); \
    \
    if (Py_TYPE(self)->tp_as_number != NULL && \
	Py_TYPE(self)->tp_as_number->NBNAME == FUNCNAME) \
	result = PyObject_CallMethodObjArgs(((cdefer_Share *)self)->runtime, \
					    name, self, other, NULL); \
    else \
	result = PyObject_CallMethodObjArgs(((cdefer_Share *)other)->runtime, \
					    name, other, self, NULL); \
    \
    Py_DECREF(name); \
    return result; \
}

SYM_BIN_OP(cdefer_Share_add, nb_add, add)
SYM_BIN_OP(cdefer_Share_mul, nb_multiply, mul)
SYM_BIN_OP(cdefer_Share_xor, nb_xor, xor)

static PyObject *cdefer_Share_sub(PyObject *self, PyObject *other) {
    PyObject* result;
    PyObject* name = PyString_FromString("sub");

    if (Py_TYPE(self)->tp_as_number != NULL &&
	Py_TYPE(self)->tp_as_number->nb_subtract == cdefer_Share_sub)
	result =  PyObject_CallMethodObjArgs(((cdefer_Share *)self)->runtime,
					     name, self, other, NULL);
    else
	result = PyObject_CallMethodObjArgs(((cdefer_Share *)other)->runtime,
					    name, self, other, NULL);

    Py_DECREF(name);
    return result;
}

static PyObject *cdefer_Share_div(PyObject *self, PyObject *other) {
    PyObject* result;
    PyObject* name = PyString_FromString("div");

    if (Py_TYPE(self)->tp_as_number != NULL &&
	Py_TYPE(self)->tp_as_number->nb_divide == cdefer_Share_div)
	result =  PyObject_CallMethodObjArgs(((cdefer_Share *)self)->runtime,
					     name, self, other, NULL);
    else
	result = PyObject_CallMethodObjArgs(((cdefer_Share *)other)->runtime,
					    name, self, other, NULL);

    Py_DECREF(name);
    return result;
}

static PyObject *cdefer_Share_floordiv(PyObject *self, PyObject *other) {
    PyObject* result;
    PyObject* name = PyString_FromString("floordiv");

    if (Py_TYPE(self)->tp_as_number != NULL &&
	Py_TYPE(self)->tp_as_number->nb_floor_divide == cdefer_Share_floordiv)
	result =  PyObject_CallMethodObjArgs(((cdefer_Share *)self)->runtime,
					     name, self, other, NULL);
    else
	result = PyObject_CallMethodObjArgs(((cdefer_Share *)other)->runtime,
					    name, self, other, NULL);

    Py_DECREF(name);
    return result;
}


#define SYM_UN_OP(FUNCNAME, SHORTNAME)	\
static PyObject *FUNCNAME(PyObject *self) { \
    PyObject* result; \
    PyObject* name = PyString_FromString(#SHORTNAME); \
    \
	result = PyObject_CallMethodObjArgs(((cdefer_Share *)self)->runtime, \
					    name, self, NULL); \
    \
    Py_DECREF(name); \
    return result; \
}

SYM_UN_OP(cdefer_Share_neg, neg)
SYM_UN_OP(cdefer_Share_invert, invert)

static PyObject *cdefer_Share_pow(PyObject *self, PyObject *other, PyObject *modulus) {
    PyObject* result, *name;

    if (Py_TYPE(self)->tp_as_number != NULL &&
	Py_TYPE(self)->tp_as_number->nb_power == cdefer_Share_pow) {
	name = PyString_FromString("pow");
	result =  PyObject_CallMethodObjArgs(((cdefer_Share *)self)->runtime,
					  name, self, other, NULL);
	Py_DECREF(name);
	return result;
    } else {
	Py_INCREF(Py_NotImplemented);
	return Py_NotImplemented;
    }
}

static PyObject *cdefer_Share_richcompare(PyObject *self,
					  PyObject *other, int op) {
    PyObject *first, *second, *tmp = NULL, *name = NULL, *one, *result;

    if (Py_TYPE(self)->tp_richcompare == cdefer_Share_richcompare) {
	first = self;
	second = other;
    } else {
	first = other;
	second = self;
    }

    switch(op) {
    case Py_EQ:
    case Py_NE:
	name = PyString_FromString("equal");
	tmp = PyObject_CallMethodObjArgs(((cdefer_Share *)first)->runtime,
					 name, first, second, NULL);
	break;

    case Py_GE:
    case Py_LT:
	name = PyString_FromString("greater_than_equal");
	tmp = PyObject_CallMethodObjArgs(((cdefer_Share *)first)->runtime,
					 name, self, other, NULL);
	break;

    case Py_LE:
    case Py_GT:
	name = PyString_FromString("greater_than_equal");
	tmp = PyObject_CallMethodObjArgs(((cdefer_Share *)first)->runtime,
					 name, other, self, NULL);
	break;
    }

    Py_XDECREF(name);

    switch(op) {
    case Py_EQ:
    case Py_GE:
    case Py_LE:
	return tmp;

    case Py_NE:
    case Py_LT:
    case Py_GT:
	if (!tmp)
	    return NULL;
	else {
	    one = PyInt_FromLong(1);
	    result = PyNumber_Subtract(one, tmp);
	    Py_DECREF(one);
	    Py_DECREF(tmp);
	    return result;
	}
    }

    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
}

static PyObject* cdefer_Share_fp(PyObject* self, void* noarg) {
    return PyInt_FromLong(((cdefer_Share*)self)->fp);
}

static PyGetSetDef cdefer_Share_getset[] = {
    {"fp", cdefer_Share_fp, 0, 0, 0},
    {0, 0, 0, 0, 0}
};

static PyNumberMethods share_as_number = {
        (binaryfunc)cdefer_Share_add,   /*nb_add*/
        (binaryfunc)cdefer_Share_sub,   /*nb_subtract*/
        (binaryfunc)cdefer_Share_mul,   /*nb_multiply*/
        (binaryfunc)cdefer_Share_div,   /*nb_divide*/
        0,                              /*nb_remainder*/
        0,                              /*nb_divmod*/
        (ternaryfunc)cdefer_Share_pow,  /*nb_power*/
/* BS       0,        */                      /*nb_negative*/
        (unaryfunc)cdefer_Share_neg,   /*nb_negative*/
        0,                              /*tp_positive*/
        0,                              /*tp_absolute*/
        0,                              /*tp_nonzero*/
/* BS        0,       */                       /*nb_invert*/
        (unaryfunc)cdefer_Share_invert, /*nb_invert*/
        0,                              /*nb_lshift*/
        0,                              /*nb_rshift*/
        0,                              /*nb_and*/
        (binaryfunc)cdefer_Share_xor,   /*nb_xor*/
        0,                              /*nb_or*/
        0,                              /*nb_coerce*/
        0,                              /*nb_int*/
        0,                              /*nb_long*/
        0,                              /*nb_float*/
        0,                              /*nb_oct*/
        0,                              /*nb_hex*/
        0,                              /* nb_inplace_add */
        0,                              /* nb_inplace_subtract */
        0,                              /* nb_inplace_multiply */
        0,                              /* nb_inplace_divide */
        0,                              /* nb_inplace_remainder */
        0,                              /* nb_inplace_power */
        0,                              /* nb_inplace_lshift */
        0,                              /* nb_inplace_rshift */
        0,                              /* nb_inplace_and */
        0,                              /* nb_inplace_xor */
        0,                              /* nb_inplace_or */
        (binaryfunc)cdefer_Share_floordiv,     /* nb_floor_divide */
        0,                              /* nb_true_divide */
        0,                              /* nb_inplace_floor_divide */
        0                               /* nb_inplace_true_divide */
};

static PyTypeObject cdefer_ShareType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "cdefer.Share",             /*tp_name*/
    sizeof(cdefer_Share),       /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)cdefer_Share_dealloc,      /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    &share_as_number,           /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE|Py_TPFLAGS_HAVE_GC|Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    0,                          /*tp_doc*/
    (traverseproc)cdefer_Share_traverse,   /*tp_traverse*/
    (inquiry)cdefer_Share_clear,           /*tp_clear*/
    (richcmpfunc)cdefer_Share_richcompare, /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    cdefer_Share_methods,       /*tp_methods*/
    cdefer_Share_members,       /*tp_members*/
    cdefer_Share_getset,         /*tp_getset*/
    &cdefer_DeferredType,       /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    (initproc)cdefer_Share_init,          /*tp_init*/
    0,                          /*tp_alloc*/
    0,                          /*tp_new*/
    PyObject_GC_Del,            /*tp_free*/
    0,                          /*tp_is_gc*/
    0,                          /*tp_bases*/
    0,                          /*tp_mro*/
    0,                          /*tp_cache*/
    0,                          /*tp_subclasses*/
    0,                          /*tp_weaklist*/
};

typedef struct {
    cdefer_Share share;
    PyObject* results;
    int missing_shares;
	long fp;
} cdefer_ShareList;

static struct PyMemberDef cdefer_ShareList_members[] = {
    {"results", T_OBJECT_EX, offsetof(cdefer_ShareList, results), 0, 0},
    {"missing_shares", T_INT, offsetof(cdefer_ShareList, missing_shares), 0, 0},
//  {"fp", T_OBJECT_EX, offsetof(cdefer_ShareList, fp), 0, 0},
    {0, 0, 0, 0, 0}
};

static void cdefer_ShareList_dealloc(PyObject *o) {
    cdefer_ShareList *self;
    self = (cdefer_ShareList *)o;
    Py_XDECREF(self->results);
    cdefer_Share_dealloc((PyObject *)&self->share);
}

static int cdefer_ShareList_traverse(PyObject *o, visitproc visit, void *arg) {
    cdefer_ShareList *self;
    self = (cdefer_ShareList *)o;
    Py_VISIT(self->results);
    return cdefer_Share_traverse((PyObject *)&self->share, visit, arg);
}

static int cdefer_ShareList_clear(PyObject *o) {
    cdefer_ShareList *self;
    self = (cdefer_ShareList *)o;
    Py_CLEAR(self->results);
    return cdefer_Share_clear((PyObject *)&self->share);
}

static int cdefer_ShareList_init(cdefer_ShareList *self, PyObject *args,
				 PyObject *kwargs)
{
    PyObject *shares, *callback_fired, *result;
    PyListObject *tmp = 0;
    cdefer_Share *first_item;
    int n_shares, threshold = 0, i, j;
    long fp_cnt = 0;
    static char *argnames[] = {"shares", "threshold", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|i", argnames,
				     &shares, &threshold))
        return -1;
	
    if (!PyList_Check(shares)) {
	tmp = (PyListObject *)PyList_New(0);
	
	if (!_PyList_Extend(tmp, shares))
	    return -1;

	shares = (PyObject *)tmp;
    }

    n_shares = PyList_Size(shares);

    if (n_shares < 0)
	return -1;

    if (n_shares < 1) {
	PyErr_SetString(PyExc_AssertionError, "Cannot create empty ShareList");
	return -1;
    }

    if (threshold < 0 || threshold > n_shares) {
	PyErr_SetString(PyExc_AssertionError, "Threshold out of range");
	return -1;
    }

    first_item = (cdefer_Share *)PyList_GET_ITEM(shares, 0);
    fp_cnt = first_item->fp;
	
/*	for (j=1; j < n_shares; j++){
		fp_cnt += ((cdefer_Share *)PyList_GET_ITEM(shares, j))->fp;
	}
	if (fp_cnt > 0){
		fp_cnt = (long)1;
	}*/

	for (j=1; j < n_shares; j++){
		fp_cnt = (long)(fp_cnt || ((cdefer_Share *)PyList_GET_ITEM(shares, j))->fp);
	}

	
    cdefer_Share__init(&(self->share), first_item->runtime,
		       first_item->field, NULL, fp_cnt);

    Py_XDECREF(self->results);
    self->results = PyList_New(n_shares);
    self->fp = fp_cnt;

    if (!self->results)
	return -1;

    if (threshold == 0)
	self->missing_shares = n_shares;
    else
	self->missing_shares = threshold;

	
    callback_fired = PyObject_GetAttrString((PyObject *)self,
					    "_callback_fired");

    if (!callback_fired)
	return -1;

    // combine loops, unlike Python code
    for (i = 0; i < n_shares; i++) {
	Py_INCREF(Py_None);
	PyList_SET_ITEM(self->results, i, Py_None);

	args = Py_BuildValue("(iO)", i, Py_True);
	Py_INCREF(Py_None);
	result = cdefer_Deferred__addCallbacks((cdefer_Deferred *)
					       PyList_GET_ITEM(shares, i),
					       callback_fired, Py_None, args,
					       Py_None, Py_None, Py_None);
	Py_DECREF(Py_None);
	Py_DECREF(args);
	Py_DECREF(result);

	if (!result)
	    return -1;
    }

    Py_DECREF(callback_fired);
    Py_XDECREF(tmp);

    return 0;
}

static PyObject *cdefer_ShareList_callback_fired(cdefer_ShareList *self,
						 PyObject *args,
						 PyObject *kwargs) {
    PyObject *result, *success, *tmp;
    int index;
    static char *argnames[] = {"result", "index", "success", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OiO", argnames,
				     &result, &index, &success))
	return NULL;

    tmp = PyList_GET_ITEM(self->results, index);
    PyList_SET_ITEM(self->results, index, Py_BuildValue("(OO)", success,
						    result));
    Py_XDECREF(tmp);
    self->missing_shares--;
  
    if (!self->share.deferred.called && self->missing_shares == 0) {
	tmp = cdefer_Deferred__startRunCallbacks(&self->share.deferred,
						 self->results);

	if (!tmp)
	    return NULL;

	Py_DECREF(tmp);
    }

    Py_INCREF(result);
    return result;
}

static PyObject* cdefer_ShareList_fp(PyObject* self, void* noarg) {
    return PyInt_FromLong(((cdefer_ShareList*)self)->fp);
}


static struct PyMethodDef cdefer_ShareList_methods[] = {
    {"_callback_fired", (PyCFunction)cdefer_ShareList_callback_fired,
     METH_VARARGS|METH_KEYWORDS, 0},
    {0, 0, 0, 0}
};



static PyGetSetDef cdefer_ShareList_getset[] = {
    {"fp", cdefer_ShareList_fp, 0, 0, 0},
    {0, 0, 0, 0, 0}
};


static PyTypeObject cdefer_ShareListType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "cdefer.ShareList",         /*tp_name*/
    sizeof(cdefer_ShareList),   /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)cdefer_ShareList_dealloc,      /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE|Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    0,                          /*tp_doc*/
    (traverseproc)cdefer_ShareList_traverse,   /*tp_traverse*/
    (inquiry)cdefer_ShareList_clear,           /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    cdefer_ShareList_methods,   /*tp_methods*/
    cdefer_ShareList_members,   /*tp_members*/
    cdefer_ShareList_getset,    /*tp_getset*/
    &cdefer_ShareType,          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    (initproc)cdefer_ShareList_init,          /*tp_init*/
    0,                          /*tp_alloc*/
    0,                          /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
    0,                          /*tp_bases*/
    0,                          /*tp_mro*/
    0,                          /*tp_cache*/
    0,                          /*tp_subclasses*/
    0,                          /*tp_weaklist*/
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initcdefer(void) {
    PyObject * m = NULL;
    PyObject * f = NULL;
    PyObject * d = NULL;
    PyObject * traceback_module = NULL;

    if (PyType_Ready(&cdefer_DeferredType) < 0) {
        return;
    }

    m = Py_InitModule3("cdefer", cdefer_methods,
                       "cdefer");

    if (!m) {
        return;
    }

    Py_INCREF(&cdefer_DeferredType);
    PyModule_AddObject(m, "Deferred", (PyObject *)&cdefer_DeferredType);

    f = PyImport_ImportModule("twisted.python.failure");
    if (!f) {
        goto Error;
    }

    failure_class = PyObject_GetAttrString(f, "Failure");
    if (!failure_class) {
        goto Error;
    }

    d = PyImport_ImportModule("twisted.internet.defer");
    if (!d) {
        goto Error;
    }
    already_called = PyObject_GetAttrString(d, "AlreadyCalledError");
    if (!already_called) {
        goto Error;
    }

    debuginfo_class = PyObject_GetAttrString(d, "DebugInfo");
    if(!debuginfo_class) {
        goto Error;
    }

    traceback_module = PyImport_ImportModule("traceback");
    if (!traceback_module) {
        goto Error;
    }

    format_stack = PyObject_GetAttrString(traceback_module, "format_stack");
    if(!format_stack) {
        goto Error;
    }

    /* VIFF */

    if (PyType_Ready(&cdefer_ShareType) < 0)
	return;

    Py_INCREF(&cdefer_ShareType);
    PyModule_AddObject(m, "Share", (PyObject *)&cdefer_ShareType);

    if (PyType_Ready(&cdefer_ShareListType) < 0)
	return;

    Py_INCREF(&cdefer_ShareListType);
    PyModule_AddObject(m, "ShareList", (PyObject *)&cdefer_ShareListType);

    return;


Error:
    Py_XDECREF(f);
    Py_XDECREF(failure_class);
    Py_XDECREF(d);
    Py_XDECREF(already_called);
    Py_XDECREF(debuginfo_class);
    Py_XDECREF(traceback_module);
    Py_XDECREF(format_stack);
}

