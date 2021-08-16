/*
 * Copyright (c) 2020-2021 Samsung Electronics Co., Ltd. All rights reserved.

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef _TVG_PICTURE_IMPL_H_
#define _TVG_PICTURE_IMPL_H_

#include <string>
#include "tvgPaint.h"
#include "tvgLoader.h"

/************************************************************************/
/* Internal Class Implementation                                        */
/************************************************************************/

struct PictureIterator : Iterator
{
    Paint* paint = nullptr;

    PictureIterator(Paint* p) : paint(p) {}

    const Paint* next() override
    {
        auto ret = paint;
        paint = nullptr;
        return ret;
    }
};


struct Picture::Impl
{
    shared_ptr<LoadModule> loader = nullptr;
    Paint* paint = nullptr;
    uint32_t *pixels = nullptr;
    Picture *picture = nullptr;
    void *rdata = nullptr;              //engine data
    float w = 0, h = 0;
    bool resizing = false;

    Impl(Picture* p) : picture(p)
    {
    }

    ~Impl()
    {
        if (paint) delete(paint);
    }

    bool dispose(RenderMethod& renderer)
    {
        bool ret = true;
        if (paint) {
            ret = paint->pImpl->dispose(renderer);
        } else if (pixels) {
            ret =  renderer.dispose(rdata);
            rdata = nullptr;
        }
        return ret;
    }

    uint32_t reload()
    {
        if (loader) {
            if (!paint) {
                if (auto p = loader->paint()) {
                    paint = p.release();
                    loader->close();
                    if (w != loader->w && h != loader->h) {
                        loader->resize(paint, w, h);
                        resizing = false;
                    }
                    if (paint) return RenderUpdateFlag::None;
                }
            }
            if (!pixels) {
                pixels = const_cast<uint32_t*>(loader->pixels());
                loader->close();
                if (pixels) return RenderUpdateFlag::Image;
            }
        }
        return RenderUpdateFlag::None;
    }

    RenderTransform resizeTransform(const RenderTransform* pTransform)
    {
        //Overriding Transformation by the desired image size
        auto sx = w / loader->w;
        auto sy = h / loader->h;
        auto scale = sx < sy ? sx : sy;

        RenderTransform tmp;
        tmp.m = {scale, 0, 0, 0, scale, 0, 0, 0, 1};

        if (!pTransform) return tmp;
        else return RenderTransform(pTransform, &tmp);
    }

    void* update(RenderMethod &renderer, const RenderTransform* pTransform, uint32_t opacity, Array<RenderData>& clips, RenderUpdateFlag pFlag)
    {
        auto flag = reload();

        if (pixels) {
            auto transform = resizeTransform(pTransform);
            rdata = renderer.prepare(*picture, rdata, &transform, opacity, clips, static_cast<RenderUpdateFlag>(pFlag | flag));
        } else if (paint) {
            if (resizing) {
                loader->resize(paint, w, h);
                resizing = false;
            }
            rdata = paint->pImpl->update(renderer, pTransform, opacity, clips, static_cast<RenderUpdateFlag>(pFlag | flag));
        }
        return rdata;
    }

    bool render(RenderMethod &renderer)
    {
        if (pixels) return renderer.renderImage(rdata);
        else if (paint) return paint->pImpl->render(renderer);
        return false;
    }

    bool viewbox(float* x, float* y, float* w, float* h) const
    {
        if (!loader) return false;
        if (x) *x = loader->vx;
        if (y) *y = loader->vy;
        if (w) *w = loader->vw;
        if (h) *h = loader->vh;
        return true;
    }

    bool size(uint32_t w, uint32_t h)
    {
        this->w = w;
        this->h = h;
        resizing = true;
        return true;
    }

    bool bounds(float* x, float* y, float* w, float* h) const
    {
        if (paint) return paint->pImpl->bounds(x, y, w, h);
        if (w) *w = this->w;
        if (h) *h = this->h;
        return true;
    }

    RenderRegion bounds(RenderMethod& renderer)
    {
        if (rdata) return renderer.region(rdata);
        if (paint) return paint->pImpl->bounds(renderer);
        return {0, 0, 0, 0};
    }

    Result load(const string& path)
    {
        if (loader) loader->close();
        bool invalid;  //Invalid Path
        loader = LoaderMgr::loader(path, &invalid);
        if (!loader) {
            if (invalid) return Result::InvalidArguments;
            return Result::NonSupport;
        }
        if (!loader->read()) return Result::Unknown;
        w = loader->w;
        h = loader->h;
        return Result::Success;
    }

    Result load(const char* data, uint32_t size, const string& mimeType, bool copy)
    {
        if (loader) loader->close();
        loader = LoaderMgr::loader(data, size, mimeType, copy);
        if (!loader) return Result::NonSupport;
        if (!loader->read()) return Result::Unknown;
        w = loader->w;
        h = loader->h;
        return Result::Success;
    }

    Result load(uint32_t* data, uint32_t w, uint32_t h, bool copy)
    {
        if (loader) loader->close();
        loader = LoaderMgr::loader(data, w, h, copy);
        if (!loader) return Result::NonSupport;
        this->w = loader->w;
        this->h = loader->h;
        return Result::Success;
    }

    Paint* duplicate()
    {
        reload();

        auto ret = Picture::gen();
        if (!ret) return nullptr;

        auto dup = ret.get()->pImpl;
        if (paint) dup->paint = paint->duplicate();

        dup->loader = loader;
        dup->pixels = pixels;
        dup->w = w;
        dup->h = h;
        dup->resizing = resizing;

        return ret.release();
    }

    Iterator* iterator()
    {
        reload();
        return new PictureIterator(paint);
    }
};

#endif //_TVG_PICTURE_IMPL_H_
