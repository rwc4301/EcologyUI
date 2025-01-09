import React from 'react'
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Users, Calendar, CreditCard, BadgeCheck } from 'lucide-react'
import App from '../App'

const data = [
  {
    name: "Jan",
    total: 15,
  },
  {
    name: "Feb",
    total: 20,
  },
  {
    name: "Mar",
    total: 18,
  },
  {
    name: "Apr",
    total: 25,
  },
  {
    name: "May",
    total: 30,
  },
  {
    name: "Jun",
    total: 28,
  },
]

export function Dashboard() {
  return (
    <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
      <Card>
        <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
          <CardTitle>Settings</CardTitle>
        </CardHeader>
        <CardContent>
          
        </CardContent>
      </Card>
      <Card className="col-span-3">
        <CardHeader>
          <CardTitle>Alpha Diversity</CardTitle>
        </CardHeader>
        <CardContent className="pl-2">
          <App />
        </CardContent>
      </Card>
    </div>
  )
}
